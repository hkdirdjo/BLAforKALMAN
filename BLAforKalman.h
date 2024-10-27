#ifndef BLAFORKALMAN_H
#define BLAFORKALMAN_H

#include <BasicLinearAlgebra.h>
    using namespace BLA;

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
class KalmanFilter
{
private:

public:
    BLA::Matrix<NUM_STATES,NUM_STATES,double> f, q, p;
    BLA::Matrix<NUM_STATES,NUM_COM,double> b;
    BLA::Matrix<NUM_STATES, 1, double > x;
    BLA::Eye<NUM_STATES, NUM_STATES,double> I;
    KalmanFilter(BLA::Matrix<NUM_STATES,NUM_STATES,double>& f, //Arguments
                 BLA::Matrix<NUM_STATES,NUM_STATES,double>& q,
                 BLA::Matrix<NUM_STATES,NUM_STATES,double>& p,
                 BLA::Matrix<NUM_STATES,NUM_COM,double>& b
                ):
                f(f), // Initializer list
                q(q),
                p(p),
                b(b)
                {
                    // Constructor; allocate resources
                    x.Fill(0);
                    p.Fill(0);
                    b.Fill(0);
                }
    void predict(const BLA::Matrix<NUM_COM,1,double>& u);
    void updateGPS(const BLA::Matrix<NUM_OBS_GPS,1,double>& zGPS, float HDOP, float VDOP);
    void updateBaro(const BLA::Matrix<NUM_OBS_BARO,1,double>& zBaro, double uncertaintyBaro);
    BLA::Matrix<NUM_STATES,1,double> state();
};

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
void KalmanFilter<NUM_STATES, NUM_OBS_GPS, NUM_OBS_BARO, NUM_COM>::predict(const BLA::Matrix<NUM_COM,1,double>& u){
    this->x = this->f * this->x + this->b * u;
    this->p = this->f * this->p * ~this->f + this->q;
}

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
void KalmanFilter<NUM_STATES, NUM_OBS_GPS, NUM_OBS_BARO, NUM_COM>::updateGPS(const BLA::Matrix<NUM_OBS_GPS,1,double>& zGPS, float HDOP, float VDOP){
    BLA::Matrix< NUM_OBS_GPS, NUM_STATES, double > h = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    BLA::Matrix<NUM_OBS_GPS, NUM_OBS_GPS, double > r = { sq(3.00*HDOP), 0.00,          0.00,
                                                         0.00,          sq(3.00*HDOP), 0.00,
                                                         0.00,          0.00,          sq(10.00*VDOP)};
    BLA::Matrix<NUM_STATES, NUM_OBS_GPS, double > k;
    BLA::Matrix<NUM_OBS_GPS, NUM_OBS_GPS, double > temp;
    BLA::Matrix<NUM_STATES, NUM_STATES, double > temp2;

  temp = h*this->p*~h + r;
  k = this->p*~h*Inverse(temp);
  this->x += k*(zGPS - h*this->x);
  temp2 = this->I - k*h;
  this->p = temp2*this->p*~temp2 + k*r*~k;
}

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
void KalmanFilter<NUM_STATES, NUM_OBS_GPS, NUM_OBS_BARO, NUM_COM>::updateBaro(const BLA::Matrix<NUM_OBS_BARO,1,double>& zBaro, double uncertaintyBaro) {
    BLA::Matrix<NUM_OBS_BARO, NUM_STATES, double > h = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    BLA::Matrix<NUM_OBS_BARO, NUM_OBS_BARO, double > r = {uncertaintyBaro};
    BLA::Matrix<NUM_STATES, NUM_OBS_BARO, double > k;
    BLA::Matrix<NUM_OBS_BARO, NUM_OBS_BARO, double > temp;
    BLA::Matrix<NUM_STATES, NUM_STATES,double > temp2;

    temp = h * this->p * ~h + r;
    k = this->p * ~h * Inverse(temp);
    this->x += k * (zBaro - h * this->x);
    temp2 = I - k * h;
    this->p = temp2 * this->p * ~temp2 + k * r * ~k;
}

template <int NUM_STATES, int NUM_OBS_GPS, int NUM_OBS_BARO, int NUM_COM>
BLA::Matrix<NUM_STATES,1,double> KalmanFilter<NUM_STATES, NUM_OBS_GPS, NUM_OBS_BARO, NUM_COM>::state() {
    return this->x;
}

#endif