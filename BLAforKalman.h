#ifndef BLAFORKALMAN_H
#define BLAFORKALMAN_H

#include <BasicLinearAlgebra.h>
    using namespace BLA;

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
class KalmanFilter
{
private:
    /* data */
    BLA::Matrix<NUM_STATES,NUM_STATES,double> f, q, p;
    BLA::Matrix<NUM_STATES,NUM_COM,double> b;
    BLA::Matrix< NUM_OBS, NUM_STATES, double > h
    BLA::Matrix<NUM_OBS, NUM_OBS, double > r
    BLA::Matrix<NUM_STATES, 1, double > x;
    BLA::Matrix<NUM_OBS_GPS, 1, double > zGPS; // posE, posN, posD
    BLA::Matrix<NUM_OBS_BARO,1, double > zBaro; // posD
public:
    KalmanFilter(BLA::Matrix<NUM_STATES,NUM_STATES,double>& f, //Arguments
                 BLA::Matrix<NUM_STATES,NUM_STATES,double>& q,
                 BLA::Matrix<NUM_STATES,NUM_STATES,double>& p,
                 BLA::Matrix<NUM_STATES,NUM_COM,double>& b,
                ):
                f(f), // Initializer list
                q(q),
                p(p),
                b(b),
                h(h),
                r(r)
                {
                    // Constructor; allocate resources
                    x = BLA::Zeros<NUM_STATES,1>();
                }
};

void predict(double dt, const BLA::Matrix<NUM_COM,1,double>& u){
    // update F matrix
    f(0,3) = deltaTime;
    f(1,4) = deltaTime;
    f(2,5) = deltaTime;

    // update B matrix
    b(0,0) = 0.5*pow(deltaTime,2);
    b(1,1) = b(0,0);
    b(2,2) = b(0,0);
    b(3,0) = deltaTime;
    b(4,1) = deltaTime;
    b(5,2) = deltaTime;

    x = f*x + b*u;
    p = f*p*~f + q;
}

void updateGPS(const BLA::Matrix<NUM_OBS_GPS,1,double>& zGPS, float HDOP, float VDOP){
    BLA::Matrix< NUM_OBS_GPS, NUM_STATES, double > h = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    BLA::Matrix<NUM_OBS_GPS, NUM_OBS_GPS, double > r = { sq(3.00*HDOP), 0.00,          0.00,
                                                         0.00,          sq(3.00*HDOP), 0.00,
                                                         0.00,          0.00,          sq(10.00*VDOP)};
    BLA::Matrix<NUM_STATES, NUM_OBS_GPS, double > k;
    BLA::Matrix<NUM_OBS_GPS, NUM_OBS_GPS, double > temp;
    BLA::Matrix<NUM_STATES, NUM_STATES, double > temp2;

    temp = h*p*~h + r;
    k = p*~h*Inverse(temp);
    x += k*(zGPS - h*x);
    temp2 = BLA::Identity<NUM_STATES,NUM_STATES>() - k*h;
    p = temp2*p*~temp2 + k*r*~k;
}

void updateBarometer(double zBaro) {
    BLA::Matrix<NUM_OBS_BARO, NUM_STATES, double > h = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    BLA::Matrix<NUM_OBS_BARO, NUM_OBS_BARO, double > r = {1.0};
    BLA::Matrix<NUM_STATES, NUM_OBS_BARO, double > k;
    BLA::Matrix<NUM_OBS_BARO, NUM_OBS_BARO, double > temp;
    BLA::Matrix<NUM_STATES, NUM_STATES,double > temp2;

    temp = h*p*~h + r;
    k = p*~h*Inverse(temp);
    x += k*(zBaro - h*x);
    temp2 = BLA::Identity<NUM_STATES,NUM_STATES>() - k*h;
    p = temp2*p*~temp2 + k*r*~k;
}

BLA::Matrix<NUM_STATES,1> state() const {
    return x;
}

#endif