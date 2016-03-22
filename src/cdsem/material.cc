/**
 * @file src/cdsem/material.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "material.h"
#include <common/constant.h>
#include <common/interpolate.h>
#include <common/spline.h>

material::material(const std::string& name, double fermi, double barrier, double density) {
	_name = name;
	_fermi = fermi;
	_barrier = barrier;
	_density = density;
}

material::material(const std::string& name, double fermi, double barrier, double bandgap, double density)
: material(name, fermi, barrier, density) {
	_bandgap = bandgap;
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
		if(P*tcs <= cit->first)
			return cit->second;
	return 0;
}

material& material::set_elastic_data(double K, const std::map<double,double>& dcs_map) {
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
	_elastic_tcs[std::log(K)] = std::log(tcs);
	for(auto cit = dcs_int_map.cbegin(); cit != dcs_int_map.cend(); cit++) {
		const double theta = cit->first;
		_elastic_dcs[std::log(K)][cumulative_dcs(theta)/tcs] = theta;
	}
	return *this;
}

material& material::set_inelastic_data(double K, const std::map<double,double>& dcs_map) {
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
	dcs_int_map[K] = 0;
	const spline cumulative_dcs = spline::linear(dcs_int_map).integrate(0);
	const double tcs = cumulative_dcs(K);
	_inelastic_tcs[std::log(K)] = std::log(tcs);
	for(auto cit = dcs_int_map.cbegin(); cit != dcs_int_map.cend(); cit++) {
		const double omega_zero = cit->first;
		_inelastic_dcs[std::log(K)][cumulative_dcs(omega_zero)/tcs] = omega_zero;
	}
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

archive::ostream& operator<<(archive::ostream& oa, const material& obj) {
	oa.put_string(obj._name);
	oa.put_float64(obj._fermi);
	oa.put_float64(obj._barrier);
	oa << obj._bandgap;
	oa.put_float64(obj._density);
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
	_put_nested_map(obj._elastic_dcs);
	_put_map(obj._inelastic_tcs);
	_put_nested_map(obj._inelastic_dcs);
	_put_nested_map(obj._ionization_tcs);
	return oa;
}

archive::istream& operator>>(archive::istream& ia, material& obj) {
	ia.get_string(obj._name);
	ia.get_float64(obj._fermi);
	ia.get_float64(obj._barrier);
	ia >> obj._bandgap;
	ia.get_float64(obj._density);
	auto _get_map = [&ia](std::map<double,double>& map) {
		map.clear();
		uint32_t n;
		ia.get_uint32(n);
		for(uint32_t i = 0; i < n; i++) {
			double x;
			ia.get_float64(x);
			ia.get_float64(map[x]);
		}
	};
	auto _get_nested_map = [&ia,&_get_map](std::map<double,std::map<double,double>>& map) {
		map.clear();
		uint32_t n;
		ia.get_uint32(n);
		for(uint32_t i = 0; i < n; i++) {
			double x;
			ia.get_float64(x);
			_get_map(map[x]);
		}
	};
	_get_map(obj._elastic_tcs);
	_get_nested_map(obj._elastic_dcs);
	_get_map(obj._inelastic_tcs);
	_get_nested_map(obj._inelastic_dcs);
	_get_nested_map(obj._ionization_tcs);
	return ia;
}