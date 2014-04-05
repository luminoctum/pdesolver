#! /usr/bin/env python2.7
from Thermodiagram import *
from InitBase import *
from pycli.ode.gridop import *
from collections import OrderedDict

class InitPrimitive(InitBase):
    def set_variables(self):
        cold_t0 = self.atmos.Tref
        #warm_t0 = find_warm_t0(self.diagram, cold_t0)
        warm_t0 = 139.69
        cold_t = voyager_T(self.diagram)
        warm_t = self.diagram.moist_dry_adiabat(warm_t0)
        cold_t = hstack([warm_t[:self.diagram.icloud[1]+1], cold_t[self.diagram.icloud[1]+1:]])
        cold_tc = self.atmos.pt_from_t(cold_t, self.diagram.ptol)
        warm_tc = self.atmos.pt_from_t(warm_t, self.diagram.ptol)
        cold_tv = self.atmos.svt_from_t(cold_t, self.diagram.ptol)
        warm_tv = self.atmos.svt_from_t(warm_t, self.diagram.ptol)

        Ptol = array([self.diagram.ptol for x in self.xaxis])
        Rdist = array([self.xaxis for y in self.yaxis]).T
        Cold_tv = array([cold_tv for x in self.xaxis])
        T = array([cold_t + (warm_t - cold_t).clip(0, 1.E6) * exp(- x**2 / (2. * self.sigma**2)) for x in self.xaxis])
        Tc = T * (self.atmos.Pref / Ptol)**(self.atmos.Rgas/self.atmos.cpvap)
        Tv = self.atmos.svt_from_t(T, Ptol)
        Eta = self.atmos.smr_from_t(T, Ptol)
        Pdry = Ptol / (1. + Eta[0] + Eta[1])
        hNH3 = Eta[0] / (1. + Eta[0]) * Ptol / self.atmos.speci[0].svp_from_t(T)
        hH2O = Eta[1] / (1. + Eta[1]) * Ptol / self.atmos.speci[1].svp_from_t(T)
        Phi = zeros((self.nx, self.ny))
        for i in range(self.ny):
            Phi[:, i] = trapz(self.atmos.grav * (Tv[:, :i + 1] - Cold_tv[:, :i + 1]) / self.atmos.Tref, self.yaxis[:i + 1], axis = 1)
        r0 = 1.E6
        Mass = array([x/r0 * exp(- self.yaxis / self.atmos.Hscale) for x in self.xaxis])
        _MassX = 1. / array([x/r0 * exp(- self.yaxis / self.atmos.Hscale) for x in self.xaxisb])
        # eliminate sigularity
        _MassX[0, :] = 0.
        _MassY = 1. / array([x/r0 * exp(- self.yaxisb / self.atmos.Hscale) for x in self.xaxis])
        T_ov_tc = array([ exp(- self.atmos.grav * self.yaxis / (self.atmos.cpvap * self.atmos.Tref)) for x in self.xaxis])
        x0 = self.xlen * (1 - 0.4)
        Absorb = zeros((self.nx, self.ny))
        for i, x in enumerate(self.xaxis):
            if x <= x0:
                Absorb[i, :] = 0
            else:
                Absorb[i, :] = self.gamma * (exp(((x - x0) / (self.xlen - x0))**4) - 1)

        self.var['radial'] = self.xaxis
        self.var['radialh'] = self.xaxisb
        self.var['zplev'] = self.yaxis
        self.var['zplevh'] = self.yaxisb
        self.var['rdist'] = Rdist
        self.var['ptol'] = Ptol
        self.var['pdry'] = Pdry
        self.var['tc'] = Tc
        self.var['phi'] = Phi
        self.var['tv0'] = Cold_tv
        self.var['mass'] = Mass
        self.var['_massX'] = _MassX
        self.var['_massY'] = _MassY
        self.var['t_ov_tc'] = T_ov_tc
        self.var['xNH3'] = Eta[0]
        self.var['xH2O'] = Eta[1]
        self.var['absorb'] = Absorb
        self.var['temp'] = T
        self.var['tempv'] = Tv
        self.var['hNH3'] = hNH3
        self.var['hH2O'] = hH2O

if __name__ == '__main__':
    #pmin, pmax = 10.E3, 17.53439E5
    pmin, pmax = 10.E3, 30E5
    nlev = 64
    nx, ny = nlev, nlev
    H2O = Water(mmr = 1.3E-2)
    NH3 = Ammonia(mmr = 4.E-4)
    gamma = 0.04
    sigma = 2.E5

    Saturn = Atmosphere(
            Pref = 1.E5,
            Tref = 134.8,
            grav = 10.44,
            xHe = 0.11,
            speci = (NH3, H2O),
            )
    diagram = Thermodiagram(
            pmin = pmin, pmax = pmax, np = nlev,
            tmin = 240, tmax = 360., nt = 400,
            atmos = Saturn,
            )
    xlen = 5.E6
    xaxisb = linspace(0, xlen, nx + 1)
    xaxis = tohalf(xaxisb)
    yaxis = diagram.zlev
    yaxisb = tohalf(diagram.zlev, ext = 'both')
    ylen = yaxisb[-1] - yaxisb[0]

    setups = OrderedDict([
            ('num_grids_in_x', nx),
            ('num_grids_in_y', ny),
            ('length_x', xlen),
            ('length_y', ylen),
            ('time_start', 0.),
            ('time_step', 5.),
            ('time_end', 180000.),
            ('steps_per_frame', 20),
            ('restart', 0),
            ('f0', 2.128E-4),
            ('absorbing_layer_coefficient', gamma),
            ('heating_width', sigma),
            ('water_mixing', H2O.mmr),
            ('ammonia_mixing', NH3.mmr),
    ])
    varlist = [
            ['time', 'di', 'time', 's'],
            ['radial', 'dk', 'radial distance', 'm'],
            ['radialh', 'dK', 'staggered radial distance', 'm'],
            ['zplev', 'dj', 'distance in log pressure coordinate', 'm'],
            ['zplevh', 'dJ', 'staggered distance in log pressure coordinate', 'm'],
            ['rdist', 'jk', 'broadcasted radial distance', 'm'],
            ['ptol', 'jk', 'broadcasted pressure', 'pa'],
            ['pdry', 'ijk', 'dry air pressure', 'pa'],
            ['uwind', 'ijk', 'zonal wind', 'm/s'],
            ['vwind', 'ijk', 'meridional wind', 'm/s'],
            ['wwind', 'ijk', 'vertical wind', 'm / s'],
            ['tc', 'ijk', 'potential temperature', 'K'],
            ['phi', 'ijk', 'geopotential height anomaly', 'm^2 / s^2'],
            ['tv0', 'jk', 'basic virtual temperature', 'K'],
            ['temp', 'ijk', 'temperature', 'K'],
            ['tempv', 'ijk', 'virtual temperature', 'K'],
            ['mass', 'jk', 'coupled mass variable = r/r0 * exp(-z / H0)', ''],
            ['_massX', 'jK', 'X staggered mass variable = (r/r0 * exp(-z / H0))**-1', ''],
            ['_massY', 'Jk', 'Y staggered mass variable = (r/r0 * exp(-z / H0))**-1', ''],
            ['t_ov_tc', 'jk', 'T over theta = exp(- gz / cpT0)', ''],
            ['xH2O', 'ijk', 'H2O mixing ratio', ''],
            ['xNH3', 'ijk', 'NH3 mixing ratio', ''],
            ['svpH2O', 'ijk', 'H2O saturation vapor pressure', ''],
            ['svpNH3', 'ijk', 'NH3 saturation vapor pressure', ''],
            ['hH2O', 'ijk', 'H2O relative humidity', ''],
            ['hNH3', 'ijk', 'NH3 relative humidity', ''],
            ['qH2O', 'ijk', 'liquid H2O mixing ratio', ''],
            ['qNH3', 'ijk', 'liquid NH3 mixing ratio', ''],
            ['absorb', 'jk', 'absorbing boundary layer coefficient', ''],
    ]

    model = InitPrimitive(setups, varlist, output = 'dynamics.nc')
    model.diagram = diagram
    model.atmos = Saturn

    model.xaxis = xaxis
    model.xaxisb = xaxisb
    model.yaxis = yaxis
    model.yaxisb = yaxisb
    model.sigma = sigma
    model.gamma = gamma
    model.initialize()
    #print model.xaxis[1] - model.xaxis[0]
    #print model.yaxis[1] - model.yaxis[0]
