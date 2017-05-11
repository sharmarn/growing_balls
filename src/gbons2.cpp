#include <cassert>
#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <cmath>
#include <chrono>
#include <map>
#include <queue>
#include <limits>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <unordered_set>
#include <CGAL/point_generators_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <libdts2/Delaunay_triangulation_s2.h>
#include <libratss/GeoCoord.h>
#include <libratss/SphericalCoord.h>
#include <libratss/ProjectS2.h>

//use the ExploredMap instead of std::map in Worker::exploreRadius
// #define EXP_RAD_USE_MAP
//#define DEBUG
//#define CGAL_DEBUG
// #define GATHER_EXP_STATS

#ifdef GATHER_EXP_STATS
	#undef EXP_RAD_USE_MAP
#endif

//data items do not have their own id, instead the priority is used as a data-id

// OUTPUT FORMAT
// #points
// then points ordered according to elimination, each with
// <id> <lat> <lon> <radius> <prio> <elimination time> <elimination partner id> 

// <id> not necessarily increasing (!)

using ratss::GeoCoord;
using ratss::SphericalCoord;

namespace { //protection namespace

using NodeID = int64_t;
using CDT = dts2::Delaunay_triangulation_with_info_s2<NodeID, void>;
using Vertex = CDT::Vertex;
using Vertex_handle = CDT::Vertex_handle;
using Vertex_circulator= CDT::Vertex_circulator;

static const double DEG_TO_RAD = 0.017453292519943295769236907684886;
static const double EARTH_RADIUS_IN_CENTIMETERS = 637279756.0856;

class randFloat
{
private:
	const uint seed = std::chrono::system_clock::now().time_since_epoch().count();
	// const uint seed=0;

	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution;
	std::function<float ()> get_rand_float;
public:
	randFloat () : generator (seed), get_rand_float (std::bind (distribution, generator)) {}
	randFloat (float a, float b) : generator (seed), distribution (a, b),
	                               get_rand_float (std::bind (distribution, generator)) {}
	float get () { return get_rand_float (); }
};

struct PointInfo {
	PointInfo() : radius(0), prio(0), id(-1), collPartnerId(-1), alive(false) {}
	PointInfo(int _id, const GeoCoord & _geo, double _prio, double _radius) :
	geo(_geo), spherical(_geo), radius(_radius), prio(_prio), id(_id),
	collPartnerId(-1), alive(true)
	{}
	PointInfo(int _id, const SphericalCoord & _spherical, double _prio, double _radius) :
	geo(_spherical), spherical(_spherical), radius(_radius), prio(_prio), id(_id),
	collPartnerId(-1), alive(true)
	{}
	PointInfo(PointInfo && other);
	PointInfo & operator=(PointInfo && other);
	GeoCoord geo;
	SphericalCoord spherical;
	double radius;
	double prio;
	Vertex_handle vh;
	int id;
	int collPartnerId;
	bool alive;
	
	void setKilled();
};

struct EOEntry {
	using CoordType = GeoCoord;
	EOEntry(int _id, const GeoCoord & _coord, double _radius, double _prio, double _time, int _collisionPartnerId) :
	id(_id), coord(_coord), radius(_radius), prio(_prio), time(_time), collisionPartnerId(_collisionPartnerId)
	{}
	int id;
	GeoCoord coord;
	double radius;
	double prio;
	double time;
	int collisionPartnerId;
};
std::ostream & operator<<(std::ostream & out, const EOEntry & eoe);

struct EO3Entry {
	using CoordType = CDT::Point;
	EO3Entry(int _id, const CDT::Point & _coord, double _radius, double _prio, double _time, int _collisionPartnerId) :
	id(_id),
	x(CDT::Geom_traits::doubleValue(_coord.x())), y(CDT::Geom_traits::doubleValue(_coord.y())), z(CDT::Geom_traits::doubleValue(_coord.z())),
	radius(_radius), prio(_prio), time(_time), collisionPartnerId(_collisionPartnerId)
	{}
	int id;
	double x;
	double y;
	double z;
	double radius;
	double prio;
	double time;
	int collisionPartnerId;
};
std::ostream & operator<<(std::ostream & out, const EO3Entry & eoe);

class EOHandlingBase {
public:
	EOHandlingBase();
	virtual ~EOHandlingBase();
	virtual std::size_t bufferSize() const = 0;
	virtual void putHeader(std::size_t pointCount) = 0;
	virtual void put(int curId, const CDT::Point & p, double time, int collisionPartnerId) = 0;
	virtual void flush()  = 0;
	virtual void reserve(std::size_t num) = 0;
	virtual void open(const std::string & outFileName) = 0;
	virtual bool is_open() = 0;
	virtual void close() = 0;
};

class EONoHandling: public EOHandlingBase {
public:
	EONoHandling() {}
	virtual ~EONoHandling() {}
	virtual std::size_t bufferSize() const { return 0; }
	virtual void putHeader(std::size_t) override {}
	virtual void put(int , const CDT::Point & , double , int ) override {}
	virtual void flush() override {}
	virtual void reserve(std::size_t) override {}
	virtual void open(const std::string &) override {}
	virtual bool is_open() { return true; }
	virtual void close() override {}
};


template<typename T_EO_TYPE>
class EOHandling: public EOHandlingBase {
public:
	using value_type = T_EO_TYPE;
	using CoordType = typename value_type::CoordType;
public:
	EOHandling(const std::vector<PointInfo> & pi);
	virtual ~EOHandling() override;
	virtual std::size_t bufferSize() const override;
	virtual void putHeader(std::size_t pointCount) override;
	virtual void put(int curId, const CDT::Point & p, double time, int collisionPartnerId) override;
	virtual void flush() override;
	virtual void reserve(std::size_t num) override;
	virtual void open(const std::string & _outFileName) override;
	virtual bool is_open() override;
	virtual void close() override;
private:
	const CoordType & coord(const CDT::Point &, int curIdx) const;
	std::string outFileName(const std::string & outFilePrefix) const;
private: //external data (keeping these in a state object would be better)
	std::ofstream out;
	const std::vector<PointInfo> & pointInfo;
	std::vector<value_type> data;
};

//The idea here is that streaming array access is way faster than dereferencing randomly allocated memory
//This is of course only true for small sets.
//Futhermore we can assume that neighbor ids are somewhat randomly distributed,
//that's why we're using a bloom filter for early termination

class ExploredVector {
protected:
	bool contains(std::vector<int>::const_iterator it, const std::vector<int>::const_iterator & end, int id) const {
		for(; it < end; ++it) {
			#ifdef GATHER_EXP_STATS
			++seekCost;
			#endif
			//early termination is faster than trying to improve branches here by reading cacheline-size many items at once
			if (*it == id) {
				return true;
			}
		}
		return false;
	}
public:
	ExploredVector() {}
	bool count(int id) const {
		return contains(d().begin(), d().end(), id);
	}
	void insert(int id) {
		m_d.emplace_back(id);
	}
	void clear() {
#ifdef GATHER_EXP_STATS
		maxFill = std::max<std::size_t>(maxFill, m_d.size());
		++clearCount;
		summedFill += m_d.size();
		
#endif
		m_d.clear();
	}
protected:
	const std::vector<int> & d() const { return m_d; }
protected:
	std::vector<int> m_d;
#ifdef GATHER_EXP_STATS
public:
	static std::size_t maxFill;
	static std::size_t clearCount;
	static std::size_t summedFill;
	static std::size_t seekCost;
#endif
};

#ifdef GATHER_EXP_STATS
std::size_t ExploredVector::maxFill = 0;
std::size_t ExploredVector::clearCount = 0;
std::size_t ExploredVector::summedFill = 0;
std::size_t ExploredVector::seekCost = 0;
#endif

class ExploredVectorWithStartOffset: public ExploredVector {
public:
	ExploredVectorWithStartOffset() {
		clear();
	}
	bool count(int id) const {
		uint8_t bf = id & 0xFF;
		return contains(d().cbegin()+m_st[bf], d().cend(), id);
	}
	void insert(int id) {
		uint8_t bf = id & 0xFF;
		ExploredVector::insert(id);
		m_st[bf] = std::min<std::size_t>(m_d.size()-1, m_st[bf]);
	}
	void clear() {
		ExploredVector::clear();
		::memset(m_st, 0xFF, 256);
	}
protected:
	uint8_t m_st[256];
};

template<typename T_BASE>
class BloomFilteredMap {
public:
	using BaseMap = T_BASE;
public:
	BloomFilteredMap() {
		clear();
	}
	bool count(int id) const {
		uint8_t bf = calcBf(id);
		if ( (m_bf[bf/8] >> (bf % 8)) & 0x1 ){
			return m_d.count(id);
		}
		return false;
	}
	void insert(int id) {
		uint8_t bf = calcBf(id);
		m_bf[bf/8] |= static_cast<uint8_t>(1) << (bf % 8);
		m_d.insert(id);
	}
	void clear() {
		m_d.clear();
		::memset(m_bf, 0, 256/8);
	}
private:
	inline uint8_t calcBf(int id) const {
		return (id >> 8) & 0xFF;
	}
private:
	BaseMap m_d;
	//bloom filter (a very stupid one at least)
	//using anything more complicated would be worse than the current implementation of a map
	//the idea here is that for a single node whose neighborhood gets explored
	//the ids of the neighbors is distributed uniformly at random and there shouldn't be too many (< 100)
	uint8_t m_bf[256/8];
};

// using ExploredMap = ExploredVectorWithStartOffset;
// using ExploredMap = BloomFilteredMap< std::set<int> >;
// using ExploredMap = BloomFilteredMap< ExploredVector >;
using ExploredMap = BloomFilteredMap< ExploredVectorWithStartOffset >;
// using ExploredMap = std::set<int>;
// using ExploredMap = ExploredVector;

struct Config {
	enum OutFormats {
		OF_NONE=0, OF_RAW=1, OF_3D=2
	};
	
	std::string inFileName;
	std::string statsOut;
	bool noOut;
	std::string outFilePrefix;
	int precision;
	int generateRandomPoints;
	uint32_t minRadius; //currently only for random points generation
	uint32_t maxRadius; //currently only for random points generation
	double maxLat;
	double minLat;
	OutFormats outFormats;
	bool simple;
	bool simpleLessRam;
	bool checkCdt;
	bool checkDuplicateProjected;
	
	
	Config();
	bool valid() const;
	void help();
	
	int parse(int argc, char ** argv);
	
	void print(std::ostream & out);
};


class Timer final {
private:
	struct timeval m_begin, m_end;
public:
	Timer();
	~Timer();
	
	void start();
	void stop();

	long usecs() const;
	double secs() const;
};


class MemUsage {
public:
	MemUsage();
	void update();
	unsigned long size,resident,share,text,lib,data,dt;
};

struct Statistics {
	double pointCreationTime;
	double cdtInsertionTime;
	double elimOrderTime;
	
	//the following measured after each function
	
	std::size_t pointCreationMemUsage;
	std::size_t cdtInsertionMemUsage;
	std::size_t elimOrderMemUsage;
};

std::ostream & operator<<(std::ostream & out, const Statistics & stats);

struct Worker
{
	///ET_UPDATE: update event; includes time and nearest neighbor
	///ET_COLLISION collision event; includes time, nearest neighbor, collision partner
	enum class EventType : int {
		
		ET_UPDATE=1, ///update event; includes time and nearest neighbor
		
		ET_COLLISION=2 ///collision event; includes time, nearest neighbor, collision partner
	};
	
	struct QueueItem {
		double time;
		int id;
		int colPartId;
		EventType evt;
		QueueItem(double _time, int _id, Worker::EventType _evt);
		QueueItem(double _time, int _id, int _colPartId, Worker::EventType _evt);
		~QueueItem();
		//event queue is a max-heap
		bool operator<(const QueueItem & other) const;
	};
	
	Worker(const Config & _cfg);

	Timer myTimer;

	Config cfg;
	Statistics stats;
	MemUsage memUsage;
	std::vector<PointInfo> pointInfo;
	ExploredMap expMap;
	
	///eventQueue holds the events that need processing
	std::priority_queue<QueueItem> eventQueue;
	
	CDT cdt;
	
	EOHandlingBase * eliminationOrder;
	
	void putHeader(std::size_t pointCount);
	
	void put(int killedIdx, int killerIdx, double time);
	void insertPoints();
	void remove(int id);
	
	const GeoCoord & geo(int id) const;
	const GeoCoord & geo(const Vertex_handle & vh) const;
	const SphericalCoord & spherical(int id) const;
	const SphericalCoord & spherical(const Vertex_handle & vh) const;
	const PointInfo & getPointInfo(const Vertex_handle & vh) const;
	PointInfo & getPointInfo(const Vertex_handle & vh);
	
	const PointInfo & getPointInfo(int id) const;
	PointInfo & getPointInfo(int id);
	const Vertex_handle & vh(int id) const;
	Vertex_handle & vh(int id);
	
	/// returns number of other points within given radius (on surface)
	///uses: geo
	int exploreRadiusNaive(double radius, int startNode) const;
	
	/// returns number of other points within given radius (on surface)
	/// uses geo, vh
	int exploreRadius(double radius, int startNodeId, std::vector< int >& results);
	
	double collisionTime(const PointInfo & pi1, const PointInfo & pi2) const;
	/// get distance to nearest neighbor (on sphere)
	double getNNDistance(int nodeId) const;
	void updateQueue(PointInfo& pi, double time);
	EventType predictCollision(const PointInfo& curPi, double time, int& resultCP, double& resultTime);
	
	int start();
	void createRandomPoints(int count, uint32_t minRadius, uint32_t maxRadius);
	void readPointsFromFile(const std::string & fileName);
	   int64_t hasIdenticalPoints();
	void prepare();
	void calcEliminationOrder();
	
	//use a n^2 log(n) algo to calculate the elimination order
	//uses n^2 storage space
	void calcEliminationOrderSimple();
	
	//use a n^2 log(n) algo to calculate the elimination order
	//uses n storage space
	void calcEliminationOrderSimpleLessRam();
};

//global function definitions

///lat and lon are in degree
double haversine(const GeoCoord & from, const GeoCoord & to) {
    double lat_arc = (from.lat- to.lat) * DEG_TO_RAD;
    double lon_arc = (from.lon - to.lon) * DEG_TO_RAD;
    double lat_h = sin(lat_arc * 0.5);
    lat_h *= lat_h;
    double lon_h = sin(lon_arc * 0.5);
    lon_h *= lon_h;
    double tmp = cos(from.lat*DEG_TO_RAD) * cos(to.lat*DEG_TO_RAD);
    return 2.0 * asin(sqrt(lat_h + tmp*lon_h));
}

double haversine(const SphericalCoord & from, const SphericalCoord & to) {
    double lat_arc = (from.theta - to.theta);
    double lon_arc = (from.phi - to.phi);
    double lat_h = ::sin(lat_arc * 0.5);
    lat_h *= lat_h;
    double lon_h = ::sin(lon_arc * 0.5);
    lon_h *= lon_h;
    double tmp = ::cos(from.theta) * ::cos(to.theta);
    return 2.0 * ::asin(::sqrt(lat_h + tmp*lon_h));
}

double distance_in_centimeters(const GeoCoord & from, const GeoCoord & to) {
	double tmp = EARTH_RADIUS_IN_CENTIMETERS*haversine(from, to);
	assert(tmp >= 0);
	return tmp;
}

double distance_in_centimeters(const SphericalCoord & from, const SphericalCoord & to) {
	double tmp = EARTH_RADIUS_IN_CENTIMETERS*haversine(from, to);
	assert(tmp >= 0);
	return tmp;
}

//
//Random point and constraint creation routines
//

// std::vector<SphericalCoord>
// getRandomPolarPoints (uint number_of_points)
// {
// 
// 	static randFloat rand_float_theta(0.002, M_PI - 0.002);
// 	static randFloat rand_float_phi(0, 2*M_PI - 0.0000000001);
// 
// 	std::vector<SphericalCoord> points;
// 	points.reserve (number_of_points);
// 
// 	for (uint i = 0; i < number_of_points; i++) {
// 		points.emplace_back (rand_float_theta.get(), rand_float_phi.get());
// 	}
// 
// 	return points;
// }

std::vector<SphericalCoord>
getRandomPolarPoints(uint number_of_points) {
	std::unordered_set<SphericalCoord> points;
	points.reserve(number_of_points);
	
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Rand = CGAL::Random_points_on_sphere_3<K::Point_3>;
	
	Rand rnd(1.0);
	ratss::ProjectS2 proj;
	
	std::cout << std::endl;
	for (; points.size() < number_of_points;) {
		K::Point_3 p = *rnd;
		++rnd;
		double theta, phi;
		proj.toSpherical(p.x(), p.y(), p.z(), theta, phi, 128);
		points.emplace(theta, phi);
		if (points.size() % 1000 == 0) {
			std::cout << '\xd' << points.size() << std::flush;
		}
	}
	std::cout << std::endl;

	return std::vector<SphericalCoord>(points.begin(), points.end());
}

PointInfo& PointInfo::operator=(PointInfo&& other) {
	geo = std::move(other.geo);
	spherical = std::move(other.spherical);
	radius = other.radius;
	prio = other.prio;
	vh = std::move(other.vh);
	id = other.id;
	collPartnerId = other.collPartnerId;
	alive = other.alive;
	return *this;
}

void PointInfo::setKilled() {
	alive = false;
	collPartnerId = -1;
	vh = Vertex_handle();
}

PointInfo::PointInfo(PointInfo&& other) :
geo(std::move(other.geo)),
spherical(std::move(other.spherical)),
radius(other.radius),
prio(other.prio),
vh(std::move(other.vh)),
id(other.id),
collPartnerId(other.collPartnerId),
alive(other.alive)
{}


std::ostream & operator<<(std::ostream & out, const EOEntry & eoe) {
	out << eoe.id << ' ' << eoe.coord.lat << ' '
		<< eoe.coord.lon << ' ' << eoe.radius << ' '
		<< eoe.prio << ' ' << eoe.time << ' ' << eoe.collisionPartnerId;
	return out;
}

std::ostream & operator<<(std::ostream & out, const EO3Entry & eoe) {
	out << eoe.id << ' '
		<< eoe.x << ' ' << eoe.y << ' ' << eoe.z << ' '
		<< eoe.radius << ' ' << eoe.prio << ' ' << eoe.time << ' ' << eoe.collisionPartnerId;
	return out;
}

EOHandlingBase::EOHandlingBase() {}
EOHandlingBase::~EOHandlingBase() {}

template<typename T_EO_TYPE>
EOHandling<T_EO_TYPE>::EOHandling(const std::vector<PointInfo> & pi) :
pointInfo(pi)
{}

template<typename T_EO_TYPE>
EOHandling<T_EO_TYPE>::~EOHandling() {
	close();
}

template<typename T_EO_TYPE>
std::size_t  EOHandling<T_EO_TYPE>::bufferSize() const {
	return data.size();
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::putHeader(std::size_t pointCount) {
	out << pointCount << '\n';
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::put(int curIdx, const CDT::Point & p, double time, int collisionPartnerId) {
// 	int id = vh->info();
	const PointInfo & pi = pointInfo.at(curIdx);
	data.emplace_back(curIdx, coord(p, curIdx), pi.radius, pi.prio, time, collisionPartnerId);
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::flush() {
	for(const value_type & eoe : data) {
		out << eoe << '\n';
	}
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::reserve(std::size_t num) {
	data.reserve(num);
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::open(const std::string & _outFileName) {
	out.open(outFileName(_outFileName));
	out.precision(std::numeric_limits<double>::digits+1);
}

template<typename T_EO_TYPE>
bool EOHandling<T_EO_TYPE>::is_open() {
	return out.is_open();
}

template<typename T_EO_TYPE>
void EOHandling<T_EO_TYPE>::close() {
	if (!out.is_open()) {
		return;
	}
	flush();
	// finally output the surviving auxiliary points
	//set collision time and coordinates of auxiliary points 
	out <<
		"0 89.999999999999972 0 9.9999999999999995e-07 2147483647 1.7976931348623157e+308 -1\n"
		"1 89.999999999999972 120 9.9999999999999995e-07 2147483647 1.7976931348623157e+308 -1\n"
		"2 89.999999999999972 -120.00000000000001 9.9999999999999995e-07 2147483647 1.7976931348623157e+308 -1\n"
		"3 -90 0 9.9999999999999995e-07 2147483647 1.7976931348623157e+308 -1\n";
	out.close();
}

template<>
const EOHandling<EOEntry>::CoordType &
EOHandling<EOEntry>::coord(const CDT::Point &, int curIdx) const {
	return pointInfo.at(curIdx).geo;
}

template<>
std::string
EOHandling<EOEntry>::outFileName(const std::string & outFilePrefix) const {
	return outFilePrefix + ".raw";
}

template<>
const EOHandling<EO3Entry>::CoordType &
EOHandling<EO3Entry>::coord(const CDT::Point & p, int) const {
	return p;
}

template<>
std::string
EOHandling<EO3Entry>::outFileName(const std::string & outFilePrefix) const {
	return outFilePrefix + ".3d";
}

Config::Config() :
noOut(false),
precision(64),
generateRandomPoints(0),
minRadius(1000),
maxRadius(6000),
maxLat(89.0),
minLat(-89.0),
outFormats(OF_NONE),
simple(false),
simpleLessRam(false),
checkCdt(true),
checkDuplicateProjected(true)
{}

void Config::help() {
	std::cout << "prg\n"
		"-gr\t<num> <min radius> <max radius>\tgenerate random points\n"
		"-o\tpath\tpath to the elimination order file\n"
		"-so\tpath to the statistics file\n"
		"-p\t<int>\tprecision to use for inital trigonometric function calculation\n"
		"--nocheck\tdon't check delaunay triangulation\n"
		"--allow-duplicates\tdon't check for duplicate projected points\n"
		"-of\t(raw|3d)\tout format, default is raw\n"
		"-noout\tdon't output anything, use this for timing\n"
		"-maxlat\t<float=89>\t ignore points above this latitude\n"
		"-minlat\t<float=89>\t ignore points below this latitude\n"
		"--simple\tuse simple n^2 log(n) algo\n"
		"--simple-less-ram\tuse simple n^3 with less memory algo\n"
		"inFile.txt"
	<< std::endl;
}

int Config::parse(int argc, char ** argv) {
	if (argc < 2) {
		return -1;
	}
	for(int i(1); i < argc; ++i) {
		std::string token(argv[i]);
		if (token == "-gr" && i+3 < argc) {
			generateRandomPoints = ::atoi(argv[i+1]);
			minRadius = ::atoi(argv[i+2]);
			maxRadius = ::atoi(argv[i+3]);
			i += 3;
		}
		else if (token == "-o" && i+1 < argc) {
			outFilePrefix = std::string(argv[i+1]);
			++i;
		}
		else if (token == "-so" && i+1 < argc) {
			statsOut = std::string(argv[i+1]);
			if (statsOut == "-") {
				statsOut.clear();
			}
			++i;
		}
		else if (token == "-p" && i+1 < argc) {
			precision = ::atoi(argv[i+1]);
			++i;
		}
		else if (token == "-of" && i+1 < argc) {
			token = std::string(argv[i+1]);
			if (token == "raw") {
				outFormats = Config::OF_RAW;
			}
			else if (token == "3d" || token == "3D") {
				outFormats = Config::OF_3D;
			}
			else {
				std::cerr << "Invalid out format: " << token << std::endl;
				return -1;
			}
			++i;
		}
		else if (token == "--allow-duplicates") {
			checkDuplicateProjected = false;
		}
		else if (token == "--nocheck") {
			checkCdt = false;
		}
		else if (token == "--simple" || token == "-simple") {
			simple = true;
		}
		else if (token == "--simple-less-ram" || token == "-simple-ess-ram") {
			simple = true;
			simpleLessRam = true;
		}
		else if (token == "-noout" || token == "--noout") {
			noOut = true;
		}
		else if (token == "-h" || token == "--help") {
			return -1;
		}
		else {
			inFileName = token;
		}
	}
	if (noOut) {
		outFormats = Config::OF_NONE;
	}
	else if(outFilePrefix.size() && outFormats == Config::OF_NONE) {
		outFormats = Config::OF_RAW;
	}
	
	return 0;
}

void Config::print(std::ostream & out) {
	out << "Mode: ";
	if (noOut) {
		out << "timing";
	}
	else {
		out << "write file";
	}
	out << '\n';
	out << "Precision: " << precision << '\n';
	out << "Maximum latitude: " << maxLat << '\n';
	out << "Minimum latitude: " << minLat << '\n';
	if (outFilePrefix.size()) {
		out << "Out file path: " << outFilePrefix << '\n';
	}
	out << "Format:";
	char sep = ' ';
	if (outFormats == OF_RAW) {
		out << sep;
		sep = ',';
		out << "raw";
	}
	if (outFormats == OF_3D) {
		out << sep;
		sep = ',';
		out << "3D";
	}
	if (sep == ' ') {
		out << sep;
		out << "none";
	}
	out << '\n';
	if (generateRandomPoints) {
		out << "#random points: " << generateRandomPoints << '\n';
		out << "min radius: " << minRadius << '\n';
		out << "max radius: " << maxRadius << '\n';
	}
	else {
		out << "In file: " << inFileName << '\n';
	}
}

bool Config::valid() const {
	return precision &&
		(generateRandomPoints || inFileName.size()) &&
		(noOut || outFilePrefix.size())
		&& (!simple || outFormats == Config::OF_RAW);
}

Timer::Timer() {
	::memset(&m_begin, 0, sizeof(struct timeval));
	::memset(&m_end, 0, sizeof(struct timeval));
}

Timer::~Timer() {}

void Timer::start() {
	gettimeofday(&m_begin, NULL);
}

void Timer::stop() {
	gettimeofday(&m_end, NULL);
}

long Timer::usecs() const {
	long mtime, seconds, useconds;
	seconds  = m_end.tv_sec  - m_begin.tv_sec;
	useconds = m_end.tv_usec - m_begin.tv_usec;
	mtime = (long)((double)((seconds) * 1000*1000 + useconds) + 0.5);
	return mtime;
}

double Timer::secs() const {
	return (double) usecs() / 1000000;
}

MemUsage::MemUsage() {
	update();
}

void MemUsage::update() {
	const char* statm_path = "/proc/self/statm";

	FILE * f = fopen(statm_path,"r");
	if(!f){
		perror(statm_path);
		return;
	}
	if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld", &size,&resident,&share,&text,&lib,&data,&dt)) {
		perror(statm_path);
	}
	fclose(f);
	//at least these four come in page-size (which is usually 4 KiB, but it would be better to retrieve the page size from sys)
	resident *= 4;
	size *= 4;
	share *= 4;
}

std::ostream & operator<<(std::ostream & out, const Statistics & stats) {
	out << "Point creation time: " << stats.pointCreationTime << '\n';
	out << "CDT insertion time: " << stats.cdtInsertionTime << '\n';
	out << "Elimination order calculation time: " << stats.elimOrderTime << '\n';
	
	out << "Point creation mem usage: " << stats.pointCreationMemUsage << '\n';
	out << "CDT insertion mem usage: " << stats.cdtInsertionMemUsage << '\n';
	out << "Elimination order calculation mem usage: " << stats.elimOrderMemUsage << '\n';
	return out;
}

Worker::QueueItem::QueueItem(double _time, int _id, Worker::EventType _evt) :
QueueItem(_time, _id, -1, _evt)
{
	assert(_evt == Worker::EventType::ET_UPDATE);
}

Worker::QueueItem::QueueItem(double _time, int _id, int _colPartId, Worker::EventType _evt) :
time(_time), id(_id), colPartId(_colPartId), evt(_evt)
{}

Worker::QueueItem::~QueueItem() {}

bool Worker::QueueItem::operator<(const Worker::QueueItem& other) const {
	//std::priority_queue is a max-heap, but we need a min-heap
	if (time == other.time) {
		if (evt == other.evt) {
			return id > other.id;
		}
		else if (evt == Worker::EventType::ET_UPDATE) { //update events have higher priority
			return false;
		}
		else {
			return true;
		}
	}
	else {
		return time > other.time;
	}
}

Worker::Worker(const Config & _cfg) : cfg(_cfg), cdt(cfg.precision), eliminationOrder(0)
{
	//the 4 initial points
	pointInfo.emplace_back(0, SphericalCoord(4.44089209850063e-16, 0), std::numeric_limits<double>::max(), 0.0000001);
	pointInfo.emplace_back(1, SphericalCoord(4.44089209850063e-16,2.0943951023932), std::numeric_limits<double>::max(), 0.0000001);
	pointInfo.emplace_back(2, SphericalCoord(4.44089209850063e-16,4.18879020478639), std::numeric_limits<double>::max(), 0.0000001);
	pointInfo.emplace_back(3, SphericalCoord(M_PI,0), std::numeric_limits<double>::max(), 0.0000001);

	if(cfg.outFilePrefix.size() && !cfg.noOut) {
		if (cfg.outFormats == Config::OF_RAW) {
			eliminationOrder = new EOHandling<EOEntry>(pointInfo);
		}
		else if (cfg.outFormats == Config::OF_3D) {
			eliminationOrder = new EOHandling<EO3Entry>(pointInfo);
		}
		else {
			throw std::runtime_error("Unsupported outfile format");
		}
	}
	else if (cfg.noOut) {
		eliminationOrder = new EONoHandling();
	}
	else {
		throw std::runtime_error("Invalid config: No out given");
	}
}
	
	
void Worker::putHeader(std::size_t pointCount) {
	eliminationOrder->putHeader(pointCount);
}

void Worker::put(int killedIdx, int killerIdx, double time) {
	assert(killedIdx == killerIdx || collisionTime(getPointInfo(killedIdx), getPointInfo(killerIdx)) == time);
	eliminationOrder->put(killedIdx, vh(killedIdx)->point(), time, killerIdx);
}

const GeoCoord& Worker::geo(int id) const {
	return getPointInfo(id).geo;
}

const GeoCoord& Worker::geo(const Vertex_handle& _vh) const {
	return getPointInfo(_vh).geo;
}

const SphericalCoord& Worker::spherical(int id) const {
	return getPointInfo(id).spherical;
}

const SphericalCoord& Worker::spherical(const Vertex_handle& _vh) const {
	return getPointInfo(_vh).spherical;
}

PointInfo& Worker::getPointInfo(int id) {
	assert(pointInfo.at(id).id == id);
	return pointInfo.at(id);
}

const PointInfo& Worker::getPointInfo(int id) const {
	return pointInfo.at(id);
}

PointInfo& Worker::getPointInfo(const Vertex_handle& _vh) {
	assert(pointInfo.at(_vh->info()).id == _vh->info());
	return getPointInfo(_vh->info());
}

const PointInfo& Worker::getPointInfo(const Vertex_handle& _vh) const {
	assert(pointInfo.at(_vh->info()).id == _vh->info());
	return getPointInfo(_vh->info());
}

Vertex_handle& Worker::vh(int id) {
	return getPointInfo(id).vh;
}

const Vertex_handle& Worker::vh(int id) const {
	return getPointInfo(id).vh;
}

void Worker::remove(int id) {
	assert(getPointInfo(id).alive);
	cdt.remove(vh(id));
	getPointInfo(id).setKilled();
}

void Worker::createRandomPoints(int count, uint32_t minRadius, uint32_t maxRadius) {
	std::cout<<"Creating random points"<<std::endl;
	std::vector<SphericalCoord> _points=getRandomPolarPoints (count);
	pointInfo.reserve(pointInfo.size() + count);
	uint32_t prio = 0;
	for(const SphericalCoord & sc :_points) {
		double radius = rand()%(maxRadius-minRadius) + minRadius;
		int id = pointInfo.size();
		pointInfo.emplace_back(id, sc, prio, radius);
		++prio;
	}
}

void Worker::readPointsFromFile(const std::string & fileName) {
	std::cout.precision(std::numeric_limits<double>::digits+1);
	std::ifstream inFile(fileName);
	
	int number_of_points;
	
	inFile>>number_of_points;
	pointInfo.reserve(pointInfo.size() + number_of_points);
	
	std::cout << "Reading "<<number_of_points<<" points from file " << fileName <<std::endl;
	double lat, lon, radius;
	int prio;
	int id = 0;
	double lat_min=1000, lat_max=-1000;
	double lon_min=1000, lon_max=-1000;
	double lat_polar_min=1000, lat_polar_max=-1000;
	double lon_polar_min=1000, lon_polar_max=-1000;
	for(int i=0; i<number_of_points; i++)
	{
		inFile>>lat>>lon>>prio>>radius;
		
		if (lat > cfg.maxLat || lat < cfg.minLat) {
			std::cout << "Removing point (" << lat << ',' << lon << ")\n";
			continue;
		}
		++id;
		
		if (lat>lat_max)
			lat_max=lat;
		if (lon>lon_max)
			lon_max=lon;
		if (lat<lat_min)
			lat_min=lat;
		if (lon<lon_min)
			lon_min=lon;
		radius*=(1024*1024);
		//std::cout.precision(17);
		//std::cout<<lat<<" "<<lon<<std::endl;
		//lat=lat+(rand()%1000)/10000000000.1;
		//lon=lon+(rand()%1000)/10000000000.1;
		// std::cout<<lat<<" "<<lon<<std::endl;
		SphericalCoord sphPoint(GeoCoord(lat,lon));
		double pol_lat=sphPoint.theta, pol_lon=sphPoint.phi;
		if (pol_lat>lat_polar_max)
			lat_polar_max=pol_lat;
		if (pol_lon>lon_polar_max)
			lon_polar_max=pol_lon;
		if (pol_lat<lat_polar_min)
			lat_polar_min=pol_lat;
		if (pol_lon<lon_polar_min)
			lon_polar_min=pol_lon;
			
		pointInfo.emplace_back(id, GeoCoord(lat, lon), prio, radius);
	}
	std::cout<<"lat was in range "<<lat_min<<" to "<<lat_max<<" and lon was in range "<<lon_min<<" to "<<lon_max<<std::endl;
	std::cout<<"Polar: lat was in range "<<lat_polar_min<<" to "<<lat_polar_max<<" and lon was in range "<<lon_polar_min<<" to "<<lon_polar_max<<std::endl;
}

int Worker::start() {
	if (cfg.generateRandomPoints)
	{
		createRandomPoints(cfg.generateRandomPoints, cfg.minRadius, cfg.maxRadius);
	}
	else
	{
		readPointsFromFile(cfg.inFileName);
	}
	
	eliminationOrder->open(cfg.outFilePrefix);
	
	if (!eliminationOrder->is_open()) {
		std::cerr << "Could not open output file" << std::endl;
		return -1;
	}
	
	prepare();
	
	if (cfg.simple) {
		if (cfg.simpleLessRam) {
			calcEliminationOrderSimpleLessRam();
		}
		else {
			calcEliminationOrderSimple();
		}
	}
	else {
		calcEliminationOrder();
	}
	
	eliminationOrder->close();
	
	if (cfg.statsOut.size()) {
		std::ofstream outStatsF(cfg.statsOut);
		if (!outStatsF.is_open()) {
			std::cerr << "Could not open stats out file " << cfg.statsOut << std::endl;
			std::cout << stats << std::endl;
		}
		else {
			outStatsF << stats << std::endl;
		}
		
	}
	else {
		std::cout << stats << std::endl;
	}
	
	return 0;
}

//stuff below is for initalization

void Worker::insertPoints() {
	std::vector< std::pair<SphericalCoord,NodeID> > points_with_id;
	points_with_id.reserve (pointInfo.size()-4);

	for(std::size_t i(4), s(pointInfo.size()); i < s; ++i) {
		const PointInfo & pi = pointInfo[i];
		if (pi.alive) {
			points_with_id.emplace_back(pi.spherical, pi.id);
		}
	}
	std::cout << "Inserting " << points_with_id.size() << " points..." << std::flush;
	myTimer.start();
	cdt.insert(points_with_id.begin(), points_with_id.end());
	myTimer.stop();
	memUsage.update();
	stats.cdtInsertionTime = myTimer.secs();
	stats.cdtInsertionMemUsage = memUsage.resident;
	std::cout << "ok" << std::endl;
	if (cfg.checkCdt) {
		std::cout << "Checking cdt..." << std::flush;
		cdt.is_valid();
		std::cout << "ok" << std::endl;
	}
}

//check if there are identical points AFTER projection
int64_t Worker::hasIdenticalPoints() {
	typedef std::pair<CDT::Point, std::size_t> PointsData;
	auto pdSmaller = [this](const PointsData & a, const PointsData & b) {
		return a.first < b.first;
	};
	auto pdEq = [this](const PointsData & a, const PointsData & b) {
		if (a.first == b.first) {
			std::cout << "Comparing identical points "
			<< a.second << " = " << b.second
			<< " -> " << this->getPointInfo(a.second).geo << " = " << this->getPointInfo(b.second).geo << std::endl;
		}
		return a.first == b.first;
	};
	std::vector<PointsData> points3D;
	points3D.reserve(pointInfo.size());
	for(std::size_t i(4), s(pointInfo.size()); i < s; ++i) {
		points3D.emplace_back( cdt.project(spherical(i)), i );
	}
	using std::sort;
	using std::unique;
	sort(points3D.begin(), points3D.end(), pdSmaller);
	auto it = unique(points3D.begin(), points3D.end(), pdEq);
	auto numIdentical = points3D.end()-it;
	return numIdentical;
}

///setup the triangulation, point indices, radii, priority indices etc.
void Worker::prepare() {
	using std::sort;
	using std::unique;
	
	//remeber: don't touch the first 4 points!
	
	//remove duplicates, they don't make any sense (collision at time 0)
	sort(pointInfo.begin()+4, pointInfo.end(), [](const PointInfo & a, const PointInfo & b) {
		if (a.geo.lat == b.geo.lat) {
			using std::abs;
			double al = abs<double>(a.geo.lon);
			double bl = abs<double>(b.geo.lon);
			if ( (a.geo.lon == b.geo.lon) || (al == bl && al == 180.0) ) {
				return false;
			}
			return a.geo.lon < b.geo.lon;
		}
		else {
			return a.geo.lat < b.geo.lat;
		}
	});
	{
		auto it = unique(pointInfo.begin()+4, pointInfo.end(), [](const PointInfo & a, const PointInfo & b) {
			if (a.geo.lat == b.geo.lat) {
				using std::abs;
				double al = abs<double>(a.geo.lon);
				double bl = abs<double>(b.geo.lon);
				bool eq = ( (a.geo.lon == b.geo.lon) || (al == bl && al == 180.0) );
				return eq;
			}
			return false;
		});
		auto numRemoved = pointInfo.end() - it;
		if (numRemoved) {
			std::cout << "Removing " << numRemoved << " identical points" << std::endl;
		}
		pointInfo.resize(it-pointInfo.begin());
	}
	//ids don't match anymore, reassign them
	for(std::size_t i(0), s(pointInfo.size()); i < s; ++i) {
		pointInfo.at(i).id = i;
	}
	
	if (cfg.checkDuplicateProjected) {
		std::cout << "Checking for points with identical projection coordinates..." << std::flush;
		int64_t numIdenticalPoints = 0;
		if (!cfg.simple && (numIdenticalPoints = hasIdenticalPoints()) ) {
			std::cout << "There are " << numIdenticalPoints << ". Existing" << std::endl;
			abort();
		}
		else {
			std::cout << "ok" << std::endl;
		}
	}

	//reserve space for our eliminationOrder
	eliminationOrder->reserve(pointInfo.size());
	
	//only the sophisticated version needs this stuff
	if (!cfg.simple) {
		cdt = CDT(cfg.precision);
		//reinit changes the vertex-handles of the first 4 points
		//also note that the first 4 points need to have an id!
		{ //assign the id of the auxiliary points
			int id = 0;
			for(auto vit(cdt.finite_vertices_begin()), vend(cdt.finite_vertices_end()); vit != vend; ++vit) {
				vit->info() = id;
				++id;
			}
		}
		insertPoints();
		
		//now set the vertex handles of the points
		for(auto vit(cdt.finite_vertices_begin()), vend(cdt.finite_vertices_end()); vit != vend; ++vit) {
			pointInfo.at(vit->info()).vh = vit;
		}
	}
}

	
//stuff below is for the calculation of elim orderes


// now check all of them for radius
int Worker::exploreRadiusNaive(double radius, int startNode) const {
	int countNearbys=0;

	const SphericalCoord & curPt = spherical(startNode);

	// iterate over all vertex handles except for auxiliary ones (NO, ALL!)
	for(auto vit(cdt.finite_vertices_begin()), vend(cdt.finite_vertices_end()); vit!=vend; ++vit) {
		const SphericalCoord & othPt = spherical(vit);
		double dist = distance_in_centimeters(curPt, othPt);
		if ((dist<=radius) && (curPt!=othPt)) {
			countNearbys++;
		}
	}
	return countNearbys;
}

int Worker::exploreRadius(double radius, int startNodeId, std::vector<int> & results) {
// idea:
// have global array for marking 'close' nodes which are within radius and one for already 'checked'
// and queue which contains nodes still to be checked 
// nodes are only put into queue, if a neighbor was checked and close
// initially, check+mark center and put all neighbors in queue

// TODO: replace map-stuff using global array with efficient reset (!) --- maybe not that important, working set for map seems small enough
//	to be cache-efficient
//For Stuttgart map is fast than local unordered_map and global unordered_map
//For Germany using the ExploredMap is faster than std::map (but not by much, generally std::map seems to be the fastest)

	results.clear();
	int countNearbys=0;

#ifdef EXP_RAD_USE_MAP
	std::map<int, bool> checked;
#else
	expMap.clear();
#endif
	std::queue<int> myQ;
	CDT::Vertex_circulator vcr, vcrEnd;

	
	myQ.push(startNodeId); // put current one into queue
	
	const SphericalCoord & centerPt = spherical(startNodeId);
	while (!myQ.empty())
	{
		int curId=myQ.front();
		myQ.pop();
		//curVhl may very well be not in the map already
		//BUT: C++ guarantees us that if it's not then it will be initalized with 0
		//See https://stackoverflow.com/questions/4523959/stdmap-default-value-for-build-in-type
	#ifdef EXP_RAD_USE_MAP
		auto & checkedRef = checked[curId];
		if (!checkedRef) {
			checkedRef=true;
	#else
		if (!expMap.count(curId)) {
			expMap.insert(curId);
	#endif
			// first check, if it is close enough
			const SphericalCoord & curPt = spherical(curId);
			double dist = distance_in_centimeters(centerPt, curPt);
			if (dist <= radius) {
				if (curPt!=centerPt) {
					countNearbys++;
					results.emplace_back(curId);
				}
				// now inspect all neighbors and put them in queue if not already checked
				vcr=cdt.incident_vertices(vh(curId));
				vcrEnd=vcr;
				if ( ! vcr.is_empty()) {
					do {
						if (cdt.is_infinite(vcr)) {
							continue;
						}
					#ifdef EXP_RAD_USE_MAP
						if (!checked[vcr->info()]) {
					#else
						if (!expMap.count(vcr->info())) {
					#endif
							myQ.push(vcr->info());
						}
					} while (++vcr != vcrEnd);
				}
			}
		}
	}
	return countNearbys;
}

double Worker::collisionTime(const PointInfo & pi1, const PointInfo & pi2) const {
	double dist = distance_in_centimeters(pi1.spherical, pi2.spherical);
	return (dist/(pi1.radius+pi2.radius));
}

// get distance to nearest neighbor (on sphere) 
double Worker::getNNDistance(int nodeId) const {
	CDT::Vertex_circulator vcr(cdt.incident_vertices(vh(nodeId)));
	CDT::Vertex_circulator vcrEnd(vcr);
	const SphericalCoord & curPt = spherical(nodeId);
	double minDist=std::numeric_limits<double>::max();
	
	if ( ! vcr.is_empty()) {
		do {
			if (cdt.is_infinite(vcr)) {
				continue;
			}
			const SphericalCoord & neiPt = spherical(vcr);
			double dist = distance_in_centimeters(curPt, neiPt);
			if (dist < minDist) {
				minDist = dist;
			}
		} while (++vcr != vcrEnd);
	}
	// inspect all neighbors, compute their distance 
	return minDist;
}

/// perform collision prediction
/// creation of update events in queue and update of NN and CP sets is NOT done in this function
Worker::EventType Worker::predictCollision(const PointInfo & curPi, double time, int & resultCP, double & resultTime) {
	assert(time >= 0);
	assert(curPi.id > 3);
	
	double dist=getNNDistance(curPi.id);
	resultTime=dist/(2*curPi.radius);
	assert(resultTime  < std::numeric_limits<double>::max());
	
	if (time<resultTime) { // goto sleep
		return EventType::ET_UPDATE;
	}
	else { // do real prediction
		std::vector<int> neighborSet;
		exploreRadius(2*dist, curPi.id, neighborSet);
		assert(neighborSet.size());
		resultTime = std::numeric_limits<double>::max();
		for(const int parId : neighborSet) {
			const PointInfo & parPi = getPointInfo(parId);
			double curTime = collisionTime(curPi, parPi); 
			if (curTime < resultTime) {
				resultTime = curTime;
				resultCP = parId;
			}
		}
		return EventType::ET_COLLISION;
	}
}

void Worker::updateQueue(PointInfo & pi, double time) {
	
	int resultCP = -1;
	double resultTime = std::numeric_limits<double>::max();

	EventType type = predictCollision(pi, time, resultCP, resultTime);
	
	if (type == EventType::ET_UPDATE) {
		pi.collPartnerId = -1;
		eventQueue.emplace(resultTime, pi.id, EventType::ET_UPDATE);
	}
	else if (type == EventType::ET_COLLISION) {
		assert(resultTime<std::numeric_limits<double>::max());
		
		pi.collPartnerId = resultCP;
		eventQueue.emplace(resultTime, pi.id, resultCP, EventType::ET_COLLISION);
	}
	else {
		assert(false);
	}
}

///the real deal, needs a valid triangulation and indices
void Worker::calcEliminationOrder () {

	putHeader(pointInfo.size());

	myTimer.start();

	// the real thing
	// invariants:
	// * for every vhl, we have predicted next collision (if any)
	// event queue contains two kinds of events (ordered according to time)
	// * collision events (stored with one of the collision partners)
	// * wake-up events (recheck for collision partners); modelled as collision event with no collision partner set
	// have an event queue which contains 

	// ADDRESS RADII, PRIOS and COLLPARTNERS via their INDEX!!!

	// now for all nodes determine nearest neighbor distance, divide that by 2*radius to get wake-up event
	// skip the first 4 aux points
	using std::next;
	for(int id(4), s(pointInfo.size()); id < s; ++id) {
		PointInfo & pi = getPointInfo(id);
		double dist = getNNDistance(id);
		eventQueue.emplace(dist/(2*pi.radius), id, EventType::ET_UPDATE);
	}
	assert(eventQueue.size() == (pointInfo.size()-4));


	while (!eventQueue.empty()) {
	
		QueueItem ce( std::move( eventQueue.top() ) );
		eventQueue.pop();
		
		PointInfo & cpi = getPointInfo(ce.id);
		
		//update event, but current point is dead -> needs no prediction
		//the 4 aux points also need no prediction
		if ((ce.evt == EventType::ET_UPDATE && !cpi.alive) || cpi.id < 4) {
			continue;
		}
		
		//this is an old collision event
		//for any type of event the collPartnerId needs to be the same
		if (ce.colPartId != cpi.collPartnerId) {
			continue;
		}

		//update event
		if (ce.evt == EventType::ET_UPDATE) {
			assert(cpi.alive);
			updateQueue(cpi, ce.time);
		}
		else {
			PointInfo & ppi = getPointInfo(ce.colPartId);
			
			if (cpi.alive && ppi.alive) {
				if (cpi.prio > ppi.prio) {
					put(ppi.id, cpi.id, ce.time);
					remove(ppi.id);
					updateQueue(cpi, ce.time);
				}
				else {
					put(cpi.id, ppi.id, ce.time);
					remove(cpi.id);
					
					//aux points do not actively collide 
					if (ppi.id < 4) {
						continue;
					}
					updateQueue(ppi, ce.time);
				}
			}
			else if (cpi.alive) { //colPart is dead
				updateQueue(cpi, ce.time);
			}
			else if (ppi.alive) { //cur is dead
				//aux points do not actively collide 
				if (ppi.id < 4) {
					continue;
				}
				updateQueue(ppi, ce.time);
			}
			else {
				; //beide tot
			}
		}
		
		if (eventQueue.size()%10000==0) {
			std::cout << '\xd' << eventQueue.size()/10000<<"0k        " << std::flush;
		}
	} //end event loop
	assert(cdt.number_of_vertices() == 4);

	std::cout << std::endl;
	
	myTimer.stop();
	stats.elimOrderTime = myTimer.secs();
	memUsage.update();
	stats.elimOrderMemUsage = memUsage.resident;
	std::cout << "Elimination took "<<myTimer.secs()<<" seconds"<<std::endl;
}

void Worker::calcEliminationOrderSimple() {
	struct QueueEntry {
		double time;
		int id1;
		int id2;
		QueueEntry(double _time, int _id1, int _id2) :
		time(_time), id1(_id1), id2(_id2)
		{}
	};
	
	putHeader(pointInfo.size());
	
	std::vector<QueueEntry> queue;
	queue.reserve((pointInfo.size()-4)*(pointInfo.size()-4));
	
	myTimer.start();
	
	std::cout << "SimpleElimination: filling the queue..." << std::flush;
	//fill the queue
	for(int id1(4), s(pointInfo.size()); id1 < s; ++id1) {
		const PointInfo & id1pi = getPointInfo(id1);
		for(int id2(id1+1); id2 < s; ++id2) {
			double ct = collisionTime(id1pi, getPointInfo(id2));
			assert(id1pi.id != id2);
			queue.emplace_back(ct, id1, id2);
		}
	}
	std::cout << "done" << std::endl;
	
	std::cout << "SimpleElimination: sorting the queue..." << std::flush;
	using std::sort;
	sort(queue.begin(), queue.end(), [](const QueueEntry & a, const QueueEntry & b) {
		if (a.time == b.time) {
			if (a.id1 == b.id1) {
				return a.id2 < b.id2;
			}
			return a.id1 < b.id1;
		}
		return a.time < b.time;
	});
	std::cout << "done" << std::endl;
	
	std::cout << "SimpleElimination: calculating order..." << std::flush;
	//queue is now sorted in ascending time order
	for(auto it(queue.begin()), end(queue.end()); it != end; ++it) {
		const QueueEntry & qe = *it;
		PointInfo & id1pi = getPointInfo(qe.id1);
		PointInfo & id2pi = getPointInfo(qe.id2);
		assert(collisionTime(id1pi, id2pi) == qe.time);
		if (!id1pi.alive || !id2pi.alive) {
			continue;
		}
		assert(id1pi.prio != id2pi.prio);
		//collision is real, kill the one with lower prio
		if (id1pi.prio < id2pi.prio) {
			put(id1pi.id, id2pi.id, qe.time);
			id1pi.alive = false;
			assert(id2pi.alive);
		}
		else {
			put(id2pi.id, id1pi.id, qe.time);
			id2pi.alive = false;
			assert(id1pi.alive);
		}
	}
	//get the last surving point
	for(int i(4), s(pointInfo.size()); i < s; ++i) {
		PointInfo & pi = getPointInfo(i);
		if (pi.alive) {
			put(pi.id, pi.id, std::numeric_limits<double>::max());
			pi.alive = false;
			break;
		}
	}
	myTimer.stop();
	memUsage.update();
	stats.elimOrderMemUsage = memUsage.resident;
	stats.elimOrderTime = myTimer.secs();
	std::cout << "done" << std::endl;
	
	std::cout << "Calculating elimination with simple algo took " << stats.elimOrderTime << " secs" << std::endl;
	
	//and we're done, now all points except for the first 4 should be dead
	#ifndef NDEBUG
	{
		int alivePoints = 0;
		for(int i(4), s(pointInfo.size()); i < s; ++i) {
			if (getPointInfo(i).alive) {
				++alivePoints;
			}
		}
		assert(!alivePoints);
	}
	#endif
}


void Worker::calcEliminationOrderSimpleLessRam() {
	struct QueueEntry {
		double time;
		int id1;
		int id2;
		QueueEntry(const QueueEntry & other) :
		time(other.time), id1(other.id1), id2(other.id2)
		{}
		QueueEntry(double _time, int _id1, int _id2) :
		time(_time), id1(_id1), id2(_id2)
		{}
		bool operator<(const QueueEntry & other) const {
			//prio queue is max-heap!
			if (time == other.time) {
				if (id1 == other.id1) {
					return id2 > other.id2;
				}
				return id1 < other.id1;
			}
			return time > other.time;
		}
	};
	
	putHeader(pointInfo.size());
	
	std::priority_queue<QueueEntry> queue;
	int remainingPoints = pointInfo.size()-4;
	
	///killedId has to be dead already!
	auto myUpdateQueue = [&queue, &remainingPoints, this](int killedId) {
		//check all points that have id as collPartnerId for their new collision partners
		//the old collision remains in the queue
		//the hope is that collisions are local an the queue will not fill up
		for(int i(4), s(pointInfo.size()); i < s; ++i) {
			PointInfo & pi = getPointInfo(i);
			if (pi.collPartnerId != killedId) {
				continue;
			}
			double ct = std::numeric_limits<double>::max();
			int collPartId = -1;
			//get next collision partner
			for(int id2(4); id2 < s; ++id2) {
				if (!getPointInfo(id2).alive || id2 == pi.id) {
					continue;
				}
				double my_ct = collisionTime(pi, getPointInfo(id2));
				if (my_ct < ct) {
					ct = my_ct;
					collPartId = id2;
				}
			}
			assert(collPartId > 3 || remainingPoints == 1);
			
			pi.collPartnerId = collPartId;
			queue.emplace(ct, pi.id, collPartId);
		}
	};
	
	myTimer.start();
	
	std::cout << "SimpleEliminationLessRam: filling the queue..." << std::flush;
	//fill the queue, find for each point the point with lowest collision time
	for(int id1(4), s(pointInfo.size()); id1 < s; ++id1) {
		PointInfo & id1pi = getPointInfo(id1);
		
		double ct = std::numeric_limits<double>::max();
		int collPartId = -1;
		//get next collision partner
		for(int id2(4); id2 < s; ++id2) {
			if (id2 == id1) {
				continue;
			}
			double my_ct = collisionTime(id1pi, getPointInfo(id2));
			if (my_ct < ct) {
				ct = my_ct;
				collPartId = id2;
			}
		}
		assert(collPartId > 3 || remainingPoints == 1);
		
		id1pi.collPartnerId = collPartId;
		queue.emplace(ct, id1, collPartId);
	}
	std::cout << "done" << std::endl;
	assert(queue.size() == (pointInfo.size()-4));
	
	std::cout << "SimpleEliminationLessRam: calculating order..." << std::flush;
	//queue is now sorted in ascending time order
	while(queue.size() && remainingPoints > 1) {
		QueueEntry qe( queue.top() );
		queue.pop();
		PointInfo & id1pi = getPointInfo(qe.id1);
		PointInfo & id2pi = getPointInfo(qe.id2);
		if (!id1pi.alive || !id2pi.alive) {
			continue;
		}
		//collision is real, kill the one with lower prio
		if (id1pi.prio < id2pi.prio) {
			put(id1pi.id, id2pi.id, qe.time);
			id1pi.alive = false;
			id1pi.collPartnerId = -1;
			--remainingPoints; //myUpdateQueue checks this
			myUpdateQueue(id1pi.id);
		}
		else {
			put(id2pi.id, id1pi.id, qe.time);
			id2pi.alive = false;
			id2pi.collPartnerId = -1;
			--remainingPoints; //myUpdateQueue checks this
			myUpdateQueue(id2pi.id);
		}
	}
	//get the last surving point
	for(int i(4), s(pointInfo.size()); i < s; ++i) {
		PointInfo & pi = getPointInfo(i);
		if (pi.alive) {
			put(pi.id, 0, std::numeric_limits<double>::max());
			pi.alive = false;
			pi.collPartnerId = -1;
			--remainingPoints;
			break;
		}
	}
	assert(remainingPoints >= 0);
	myTimer.stop();
	memUsage.update();
	stats.elimOrderMemUsage = memUsage.resident;
	stats.elimOrderTime = myTimer.secs();
	std::cout << "done" << std::endl;
	
	std::cout << "Calculating elimination with simple-less-ram algo took " << stats.elimOrderTime << " secs" << std::endl;
	
	//and we're done, now all points except for the first 4 should be dead, and the last suriving one
	#ifndef NDEBUG
	{
		int alivePoints = 0;
		for(int i(4), s(pointInfo.size()); i < s; ++i) {
			if (getPointInfo(i).alive) {
				++alivePoints;
			}
		}
		assert(alivePoints == 0);
	}
	#endif
}

//
// main
//

} //end protection namespace

int
main(int argc, char* argv[]) {
	Config cfg;
	if (cfg.parse(argc, argv) < 0) {
		cfg.help();
		return -1;
	}
	if (!cfg.valid()) {
		std::cout << "Invalid commandline options" << std::endl;
		cfg.help();
		return -1;
	}
	cfg.print(std::cout);
	
	Worker w(cfg);
	if (w.start() < 0) {
		return -1;
	}

	#ifdef GATHER_EXP_STATS
	std::cout << "Explored-map max fill: " << ExploredVector::maxFill << std::endl;
	std::cout << "Explored-map mean fill: " << (ExploredVector::summedFill/ExploredVector::clearCount) + 1 << std::endl;
	std::cout << "Explored-map seek cost: " << ExploredVector::seekCost << std::endl;
	#endif
	
	return 0;
}

