/*
 * uFVM-CPP:
 *      Parallel unstructured 3D fluid flow CFD package 2020
 *
 * Developed by:
 *      Fadl Moukalled
 *      Mhamad Mahdi Alloush
 *      Adam Fares
 *
 * File:
 *      realm.h
*/

#ifndef realm_h
#define realm_h

#include <vector>
#include <unordered_map>
#include <ostream>
#include <cmath>

namespace ufvm
{

typedef double scalar;
typedef int label;

class Point
{
private:
    double x_,y_,z_;

public:
    Point() : x_(0), y_(0), z_(0) {}
    Point(double x, double y, double z) :x_(x), y_(y), z_(z) {}
    double& x() { return x_; }
    double& y() { return y_; }
    double& z() { return z_; }

    friend std::ostream& operator<< (std::ostream& os, Point& p) {
        os << "(" << p.x() << "," << p.y() << "," << p.z() << ")";
        return os;
    }

    friend Point operator* (scalar scale, Point& p){
        return Point(scale * p.x(), scale * p.y(), scale * p.z());
    }

    friend Point operator* (Point& p, scalar scale){
        return Point(scale * p.x(), scale * p.y(), scale * p.z());
    }

    operator std::string() const {
        return std::to_string(x_) + "," + std::to_string(y_) + "," + std::to_string(z_);
    }

    scalar magnitude() const {
        return sqrt(x_ * x_ + y_ * y_ + z_*z_);
    }
};

struct Velocity
{
private:
    double _mag;
    double _dir;

public:
    Velocity(double mag, double dir) : _mag(mag), _dir(dir){}
    scalar magnitude() const {
        return _mag;
    }
    scalar direction() const {
        return _dir;
    }
};

typedef std::vector<Point> coordinates;
typedef std::vector<Velocity> velocities;
class realm
{

private:

    label globalIter_ = 0;

    label iter_ = 0;

    scalar time_ = 0;

    label timeStepCount_ = 0;

    bool plotRes_ = false;

    bool printScales_ = false;

    scalar startTime_ = 0;

    // Input Options

    scalar initial_pollution_rate; // g/s
    label passquill_stability;
    Point stack_outlet;
    coordinates receiver_points;
    std::vector<scalar> receiver_measurements;
    Point gps_uncertainty;
    scalar sensor_uncertainty;
    velocities wind_velocity;
    scalar ppm_converter_factor;
    scalar flow_rate;

    std::vector<std::pair<std::string, std::string> > reverse_model_values;
    std::vector<std::pair<std::string, std::string> > forward_model_values;
    Point distance(scalar lat1, scalar lon1, scalar lat2, scalar lon2);


public:

    // Constructors

    realm();


    // Execution

    void run();

    

    // Operations

    int signum(int num);

    void computeAlpha_x(scalar& phi_x, Point U);

    Point rotatePointInv(Point& v, Point& stack, scalar angle);

    scalar compute_Concentration_old(Point cell, scalar Q, Point U, int& stabilityClass, Point& stack_outlet);

    scalar compute_Concentration(const Point cell, scalar Q, Velocity U, int stabilityClass, Point& stack_outlet);

    scalar compute_amountQ(Point& cell, scalar C_meas, Velocity& U, int stabilityClass, Point& stack_outlet, Point GPS_unc, scalar Sensor_unc);

    void computeSigmaYSigmaZ(scalar dx, int stability, scalar& sigma_y, scalar& sigma_z, scalar& phi);

    void addrecevier_measurement(scalar val) { receiver_measurements.emplace_back(val);}
    void addreceiver_point(Point& point) { receiver_points.emplace_back(point); }
    void addwind_velocity(Velocity& vel) { wind_velocity.emplace_back(vel); }

    std::vector<std::pair<std::string, std::string> > getFlowRates() { return reverse_model_values; }
    std::vector<std::pair<std::string, std::string> > getConcentrations() { return forward_model_values; }
    void setStackOutlet(Point& point) { stack_outlet = point; }
    void setFlowRate(scalar val) { flow_rate = val; }
    //void computeSigmaYSigmaZ_phi(scalar dx, int stability, scalar& sigma_y, scalar& sigma_z, scalar& phi);

};

} // accel


#endif // realm_h
