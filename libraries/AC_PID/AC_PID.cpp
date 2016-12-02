// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-

/// @file	AC_PID.cpp
/// @brief	Generic PID algorithm

#include <AP_Math/AP_Math.h>
#include "AC_PID.h"

const AP_Param::GroupInfo AC_PID::var_info[] PROGMEM = {
    // @Param: P
    // @DisplayName: PID Proportional Gain
    // @Description: P Gain which produces an output value that is proportional to the current error value
    AP_GROUPINFO("P",    0, AC_PID, _kp, 0),

    // @Param: I
    // @DisplayName: PID Integral Gain
    // @Description: I Gain which produces an output that is proportional to both the magnitude and the duration of the error
    AP_GROUPINFO("I",    1, AC_PID, _ki, 0),

    // @Param: D
    // @DisplayName: PID Derivative Gain
    // @Description: D Gain which produces an output that is proportional to the rate of change of the error
    AP_GROUPINFO("D",    2, AC_PID, _kd, 0),

    // 3 was for uint16 IMAX
    // 4 is used by TradHeli for FF

    // @Param: IMAX
    // @DisplayName: PID Integral Maximum
    // @Description: The maximum/minimum value that the I term can output
    AP_GROUPINFO("IMAX", 5, AC_PID, _imax, 0),

    // @Param: FILT_HZ
    // @DisplayName: PID Input filter frequency in Hz
    // @Description: Input filter frequency in Hz
    // @Unit: Hz
    AP_GROUPINFO("FILT_HZ", 6, AC_PID, _filt_hz, AC_PID_FILT_HZ_DEFAULT),

    AP_GROUPEND
};

// Constructor
AC_PID::AC_PID(float initial_p, float initial_i, float initial_d, float initial_imax, float initial_filt_hz, float dt) :
    _dt(dt),
    _integrator(0.0f),
    _input(0.0f),
    _input1(0.0f),
    _input2(0.0f),
    _signal(0.0f),
    _signal1(0.0f),
    _signal2(0.0f),
    _ntchsig(0.0f),
    _ntchsig1(0.0f),
    _ntchsig2(0.0f),
    _derivative(0.0f)
{
    // load parameter values from eeprom
    AP_Param::setup_object_defaults(this, var_info);

    _kp = initial_p;
    _ki = initial_i;
    _kd = initial_d;
    _imax = fabsf(initial_imax);
    filt_hz(initial_filt_hz);

    // reset input filter to first value received
    _flags._reset_filter = true;

    memset(&_pid_info, 0, sizeof(_pid_info));
}

// set_dt - set time step in seconds
void AC_PID::set_dt(float dt)
{
    // set dt and calculate the input filter alpha
    _dt = dt;
}

// filt_hz - set input filter hz
void AC_PID::filt_hz(float hz)
{
    _filt_hz.set(fabsf(hz));

    // sanity check _filt_hz
    _filt_hz = max(_filt_hz, AC_PID_FILT_HZ_MIN);
}

// set_input_filter_all - set input to PID controller
//  input is filtered before the PID controllers are run
//  this should be called before any other calls to get_p, get_i or get_d
void AC_PID::set_input_filter_all(float input)
{
    // don't process inf or NaN
    if (!isfinite(input)) {
        return;
    }

    // reset input filter to value received
    if (_flags._reset_filter) {
        _flags._reset_filter = false;
        _input = input;
        _input1 = input;
        _input2 = input;
        _signal = input;
        _signal1 = input;
        _signal2 = input;
        _ntchsig = input;
        _ntchsig1 = input;
        _ntchsig2 = input;
        _derivative = 0.0f;
    }

// Original first order filter
    // update filter and calculate derivative
//    float input_filt_change = get_filt_alpha() * (input - _input);
//    _input = _input + input_filt_change;

    _ntchsig=input;
// 2nd order notch filter
    float notch_hz = 6.2f;
    if (notch_hz > 0.0f && _filt_hz < 10.0f) {
      float notch_Q = 0.8f;
      float notch_c = 1/tanf(M_PI_F*_dt*notch_hz);
      float n0 = notch_Q*(notch_c*notch_c+1.0f);
      float n1 =-2.0f*notch_Q*(notch_c*notch_c-1.0f);
      float n2 = n0;
      float d0 = notch_c*notch_c*notch_Q+notch_c+notch_Q;
      float d1 = -2.0f*notch_Q*(notch_c*notch_c-1.0f);
      float d2 = notch_c*notch_c*notch_Q-notch_c+notch_Q;
      _signal = (n0*_ntchsig+n1*_ntchsig1+n2*_ntchsig2-d1*_signal1-d2*_signal2)/d0;
    } else {
      _signal = _ntchsig;
      _signal1 = _ntchsig1;
      _signal2 = _ntchsig2;
    }


// 2nd order butterworth filter
//    _signal = input;
    float ita = 1/tanf(M_PI_F*_dt*_filt_hz);
    float b0 = 1.0f/(1.0f+safe_sqrt(2.0f)*ita+ita*ita);
    float b1 = 2.0f*b0;
    float b2 = b0;
    float a1 = 2.0f*(ita*ita-1.0f)*b0;
    float a2 = -1.0f*(1.0f-safe_sqrt(2.0f)*ita+ita*ita)*b0;
    _input = b0*_signal+b1*_signal1+b2*_signal2+a1*_input1+a2*_input2;

// 1st order butterworth filter
//    _signal = input;
//    float ita = 1/tanf(M_PI_F*_dt*_filt_hz);
//    float b0 = 1.0f/(ita+1.0f);
//    float b1 = b0;
//    float a1 = (1.0f-ita)/(ita+1.0f);
//    _input = b0*_signal+b1*_signal1+a1*_input1;


    float input_filt_change = _input - _input1;

    if (_dt > 0.0f) {
        _derivative = input_filt_change / _dt;
    }
// Update n-1, and n-2 values of filter and signal for next iteration
    _ntchsig2 = _ntchsig1;
    _ntchsig1 = _ntchsig;
    _signal2 = _signal1;
    _signal1 = _signal;
    _input2 = _input1;
    _input1 = _input;

}

// set_input_filter_d - set input to PID controller
//  only input to the D portion of the controller is filtered
//  this should be called before any other calls to get_p, get_i or get_d
void AC_PID::set_input_filter_d(float input)
{
    // don't process inf or NaN
    if (!isfinite(input)) {
        return;
    }

    // reset input filter to value received
    if (_flags._reset_filter) {
        _flags._reset_filter = false;
        _derivative = 0.0f;
    }

    // update filter and calculate derivative
    if (_dt > 0.0f) {
        float derivative = (input - _input) / _dt;
        _derivative = _derivative + get_filt_alpha() * (derivative-_derivative);
    }

    _input = input;
}

float AC_PID::get_p()
{
    _pid_info.P = (_input * _kp);
    return _pid_info.P;
}

float AC_PID::get_i()
{
    if(!is_zero(_ki) && !is_zero(_dt)) {
        _integrator += ((float)_input * _ki) * _dt;
        if (_integrator < -_imax) {
            _integrator = -_imax;
        } else if (_integrator > _imax) {
            _integrator = _imax;
        }
        _pid_info.I = _integrator;
        return _integrator;
    }
    return 0;
}

float AC_PID::get_d()
{
    // derivative component
    _pid_info.D = (_kd * _derivative);
    return _pid_info.D;
}

float AC_PID::get_pi()
{
    return get_p() + get_i();
}

float AC_PID::get_pid()
{
    return get_p() + get_i() + get_d();
}

void AC_PID::reset_I()
{
    _integrator = 0;
}

void AC_PID::load_gains()
{
    _kp.load();
    _ki.load();
    _kd.load();
    _imax.load();
    _imax = fabsf(_imax);
    _filt_hz.load();
}

// save_gains - save gains to eeprom
void AC_PID::save_gains()
{
    _kp.save();
    _ki.save();
    _kd.save();
    _imax.save();
    _filt_hz.save();
}

/// Overload the function call operator to permit easy initialisation
void AC_PID::operator() (float p, float i, float d, float imaxval, float input_filt_hz, float dt)
{
    _kp = p;
    _ki = i;
    _kd = d;
    _imax = fabsf(imaxval);
    _filt_hz = input_filt_hz;
    _dt = dt;
}

// calc_filt_alpha - recalculate the input filter alpha
float AC_PID::get_filt_alpha() const
{
    if (is_zero(_filt_hz)) {
        return 1.0f;
    }

    // calculate alpha
    float rc = 1/(M_2PI_F*_filt_hz);
    return _dt / (_dt + rc);
}
