/**
 * @file src/cdsem/material.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "material.hh"
#include <algorithm>
#include <functional>
#include "constant.hh"
#include "interpolate.hh"
#include "spline.hh"

material::material(const std::string& name, double fermi, double barrier, double phonon_loss, double density) {
    _name = name;
    _fermi = fermi;
    _barrier = barrier;
    _phonon_loss = phonon_loss;
    _density = density;
}

material::material(const std::string& name, double fermi, double barrier, double phonon_loss, double band_gap, double density)
: material(name, fermi, barrier, phonon_loss, density) {
    _band_gap = band_gap;
}

double material::ionization_energy(double K, double P) const {
    std::map<double,double> ionization_map;
    double tcs = 0;
    for(auto cit = _ionization_tcs.cbegin(); cit != _ionization_tcs.cend(); cit++)
        if(K > cit->first) {
            tcs += std::exp(interpolate(cit->second, std::log(K)));
            ionization_map[tcs] = cit->first;
        }
    for(auto cit = ionization_map.cbegin(); cit != ionization_map.cend(); cit++)
        if(P*tcs < cit->first)
            return cit->second;
    return 0;
}

double material::outer_shell_ionization_energy(double omega0) const {
    for(const double binding_energy : _osi_energies)
        if((binding_energy < 100.0*constant::ec) && (omega0 > binding_energy))
            return binding_energy;
    return -1;
}

material& material::set_elastic_data(double K, const std::map<double,double>& dcs_map) {
    const double log_K = std::log(K);
    std::map<double,double> dcs_int_map;
    for(auto cit = dcs_map.cbegin(); cit != dcs_map.cend(); cit++) {
        const double theta = cit->first;
        const double dcs = cit->second;
        if((theta > 0) && (theta < constant::pi) && (dcs > 0))
            dcs_int_map[theta] = dcs*2.0*constant::pi*std::sin(theta);
    }
    if(dcs_int_map.empty())
        return *this;
    dcs_int_map[0] = 0;
    dcs_int_map[constant::pi] = 0;
    const spline cumulative_dcs = spline::linear(dcs_int_map).integrate(0);
    const double tcs = cumulative_dcs(constant::pi);
    _elastic_tcs[log_K] = std::log(tcs);
    std::map<double,double> icdf_map;
    for(auto cit = dcs_int_map.cbegin(); cit != dcs_int_map.cend(); cit++) {
        const double theta = cit->first;
        const double P = cumulative_dcs(theta)/tcs;
        if(icdf_map.count(P) == 0)
            icdf_map.insert(std::make_pair(P, theta));
    }
    _elastic_icdf.insert(std::make_pair(log_K, icdf_map));
    return *this;
}

material& material::set_inelastic_data(double K, const std::map<double,double>& dcs_map) {
    const double log_K = std::log(K);
    std::map<double,double> dcs_int_map;
    for(auto cit = dcs_map.cbegin(); cit != dcs_map.cend(); cit++) {
        const double omega_zero = cit->first;
        const double dcs = cit->second;
        if((omega_zero > 0) && (omega_zero < K) && (dcs > 0))
            dcs_int_map[omega_zero] = dcs;
    }
    if(dcs_int_map.empty())
        return *this;
    dcs_int_map[0] = 0;
   	dcs_int_map[std::min(K, dcs_map.crbegin()->first)] = 0;
    double omega_max = std::min(K, dcs_map.crbegin()->first);
    const spline cumulative_dcs = spline::linear(dcs_int_map).integrate(0);
    const double tcs = cumulative_dcs(omega_max);
    _inelastic_tcs[log_K] = std::log(tcs);
    std::map<double,double> icdf_map;
    for(auto cit = dcs_int_map.cbegin(); cit != dcs_int_map.cend(); cit++) {
        const double omega_zero = cit->first;
        const double P = cumulative_dcs(omega_zero)/tcs;
        if(icdf_map.count(P) == 0)
            icdf_map.insert(std::make_pair(P, std::log(omega_zero)));
    }
    _inelastic_icdf.insert(std::make_pair(log_K, icdf_map));
    return *this;
}

material& material::set_ionization_data(double B, const std::map<double,double>& tcs_map) {
    std::map<double,double> loglog_tcs_map;
    for(auto cit = tcs_map.cbegin(); cit != tcs_map.cend(); cit++) {
        const double K = cit->first;
        const double tcs = cit->second;
        if((K > B) && (tcs > 0))
            loglog_tcs_map[std::log(K)] = std::log(tcs);
    }
    if(loglog_tcs_map.empty())
        return *this;
    _ionization_tcs[B] = loglog_tcs_map;
    return *this;
}

material& material::set_outer_shell_ionization_data(const std::vector<double>& osi_vector) {
    _osi_energies = osi_vector;
    std::sort(_osi_energies.begin(), _osi_energies.end(), std::greater<double>());
    return *this;
}

archive::ostream& operator<<(archive::ostream& oa, const material& obj) {
    oa.put_string(obj._name);
    oa.put_float64(obj._fermi);
    oa.put_float64(obj._barrier);
    oa << obj._band_gap;
    oa.put_float64(obj._phonon_loss);
    oa.put_float64(obj._density);
    auto _put_vector = [&oa](const std::vector<double>& vector) {
        oa.put_uint32(vector.size());
        for(auto cit = vector.cbegin(); cit != vector.cend(); cit++) {
            oa.put_float64(*cit);
        }
    };
    auto _put_map = [&oa](const std::map<double,double>& map) {
        oa.put_uint32(map.size());
        for(auto cit = map.cbegin(); cit != map.cend(); cit++) {
            oa.put_float64(cit->first);
            oa.put_float64(cit->second);
        }
    };
    auto _put_nested_map = [&oa,&_put_map](const std::map<double,std::map<double,double>>& map) {
        oa.put_uint32(map.size());
        for(auto cit = map.cbegin(); cit != map.cend(); cit++) {
            oa.put_float64(cit->first);
            _put_map(cit->second);
        }
    };
    _put_map(obj._elastic_tcs);
    _put_nested_map(obj._elastic_icdf);
    _put_map(obj._inelastic_tcs);
    _put_nested_map(obj._inelastic_icdf);
    _put_nested_map(obj._ionization_tcs);
    _put_vector(obj._osi_energies);
    return oa;
}

archive::istream& operator>>(archive::istream& ia, material& obj) {
    ia.get_string(obj._name);
    ia.get_float64(obj._fermi);
    ia.get_float64(obj._barrier);
    ia >> obj._band_gap;
    ia.get_float64(obj._phonon_loss);
    ia.get_float64(obj._density);
    auto _get_vector = [&ia](std::vector<double>& vector) {
        vector.clear();
        uint32_t n;
        ia.get_uint32(n);
        vector.resize(n);
        for(uint32_t i = 0; i < n; i++) {
            ia.get_float64(vector[i]);
        }
    };
    auto _get_map = [&ia](std::map<double,double>& map) {
        map.clear();
        uint32_t n;
        ia.get_uint32(n);
        for(uint32_t i = 0; i < n; i++) {
            double x, y;
            ia.get_float64(x);
            ia.get_float64(y);
            map.insert(std::make_pair(x, y));
        }
    };
    auto _get_nested_map = [&ia,&_get_map](std::map<double,std::map<double,double>>& map) {
        map.clear();
        uint32_t n;
        ia.get_uint32(n);
        for(uint32_t i = 0; i < n; i++) {
            double x;
            ia.get_float64(x);
            std::map<double,double> nested_map;
            _get_map(nested_map);
            map.insert(std::make_pair(x, nested_map));
        }
    };
    _get_map(obj._elastic_tcs);
    _get_nested_map(obj._elastic_icdf);
    _get_map(obj._inelastic_tcs);
    _get_nested_map(obj._inelastic_icdf);
    _get_nested_map(obj._ionization_tcs);
    _get_vector(obj._osi_energies);

    return ia;
}
