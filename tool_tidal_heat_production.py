import numpy as NP
import numpy.linalg as NPL
import scipy
from scipy.interpolate import CubicSpline as CSINTP
import scipy.special as SSP
import sympy as SY
import sympy.matrices as SYM
import tool_layerallocation as LA

class Layer(object):
    def __init__( self, inner_radius, outer_radius, 
                 Density, Viscosity, ShearWave_Velocity, PressureWave_Velocity, 
                 Heat_Capacity, radiogenic_fraction):
        self.Ri = inner_radius
        self.Ro = outer_radius
        self.thickness = outer_radius-inner_radius
        assert( self.thickness > 0.0 )
        self.Si = 4.0 * NP.pi * self.Ri**2.
        self.So = 4.0 * NP.pi * self.Ro**2.
        self.Volume = 4.0/3.0 * NP.pi * (self.Ro**3. - self.Ri**3.) 
        self.Rho = Density
        self.Eta = Viscosity
        self.Cp = Heat_Capacity
        self.Vs = ShearWave_Velocity
        self.Vp = PressureWave_Velocity
        self.radiogenic_fraction = radiogenic_fraction

def temporal_evolution_of_tidal_heating(EM_Distance, model, 
                                        Crust_Layer, Mantle_Layer, LVZ_Layer, Core_Layer):
    if (EM_Distance<17.9) or (EM_Distance>60.0):
        return 0.0, 0.0
    
    Angular_Velocity = 2 * NP.pi / ((EM_Distance/60.3)**(1.5)*29.530*24*60*60)## a~35 Re at 4.0 Ga, ~60 Re at present. According to Kepler'law, 
    Eccentricity = 0.0549
    G_Constant = 6.67e-11
    Alpha = 0.5 
    Rs = Crust_Layer.Ro
    Eta = [Core_Layer.Eta, LVZ_Layer.Eta, Mantle_Layer.Eta, Crust_Layer.Eta]
    Vs = [Core_Layer.Vs, LVZ_Layer.Vs, Mantle_Layer.Vs, Crust_Layer.Vs] 
    Vp = [Core_Layer.Vp, LVZ_Layer.Vp, Mantle_Layer.Vp, Crust_Layer.Vp] 
    Rho = [Core_Layer.Rho, LVZ_Layer.Rho, Mantle_Layer.Rho, Crust_Layer.Rho] 
    Mu = [x*y**2 for (x,y) in zip(Rho, Vs)]
    Kappa = [x*y**2 - 4*z/3 for (x,y,z) in zip(Rho, Vp, Mu)]
    Layer_Number = 101 ##Innermost layer: #0; Outermost layer: #100. We assume physical properties below are constant in any layer
    nr = [1, 41, 96, 101] ## Changed
    rr = [Core_Layer.Ro, LVZ_Layer.Ro, Mantle_Layer.Ro, Crust_Layer.Ro]
    Radius = LA.Generation_radius(nr, rr)
    RADIUS = NP.array(Radius, dtype=complex).reshape(-1,1) ## RADIUS[0] is the top radius of Layer #0, i.e., the core, so do other layers. 
    GRAVITY = NP.zeros((Layer_Number,1), dtype=complex) ## GRAVITY is the gravitational acceleration of Layer #0, so do other layers.
    DENSITY = NP.zeros((Layer_Number,1), dtype=complex) ## DENSITY is the density of Layer #0, so do other layers.
    SHEAR_MODULUS = NP.zeros((Layer_Number,1), dtype=complex) ## complex number value in the viscoelastic problem
    LAME_CONSTANT = NP.zeros((Layer_Number,1), dtype=complex) ## complex number value in the viscoelastic problem
    BULK_MODULUS = NP.zeros((Layer_Number,1), dtype=complex) ## still real number value in the viscoelastic problem
    mu_tmp = 0 ## used to calculate real/imaginary-number-value shear modulus
    sm_real = 0 ## real part of complex shear modulus
    sm_img = 0 ## imaginary part of complex shear modulus
    if model=='Maxwell':
        for layer, radius in enumerate(RADIUS):
            if radius.real <= Radius[0]:
                DENSITY[layer] = Rho[0]
                GRAVITY[layer] = 4/3*NP.pi*G_Constant * Rho[0] * radius
                mu_tmp = Mu[0]
                eta_tmp = Eta[0]
                BULK_MODULUS[layer] = Kappa[0]
                sm_real = mu_tmp * Angular_Velocity**2 * eta_tmp**2 / (mu_tmp**2 + Angular_Velocity**2 * eta_tmp**2)
                sm_img = mu_tmp**2 * Angular_Velocity * eta_tmp / (mu_tmp**2 + Angular_Velocity**2 * eta_tmp**2)
                SHEAR_MODULUS[layer] = sm_real + sm_img * 1j
                LAME_CONSTANT[layer] = BULK_MODULUS[layer] - 2*SHEAR_MODULUS[layer]/3
            else:
                layer_material = LA.Find_MaterialLayer(radius.real, rr)
                DENSITY[layer] = Rho[layer_material]
                GRAVITY[layer] = 4/3*NP.pi*G_Constant*DENSITY[layer]*radius +  GRAVITY[layer-1]*(RADIUS[layer-1]/radius)**2 - 4/3*NP.pi*G_Constant*DENSITY[layer]*RADIUS[layer-1]**3/radius**2 
                mu_tmp = Mu[layer_material]
                eta_tmp = Eta[layer_material]
                BULK_MODULUS[layer] = Kappa[layer_material]
                sm_real = mu_tmp * Angular_Velocity**2 * eta_tmp**2 / (mu_tmp**2 + Angular_Velocity**2 * eta_tmp**2)
                sm_img = mu_tmp**2 * Angular_Velocity * eta_tmp / (mu_tmp**2 + Angular_Velocity**2 * eta_tmp**2)
                SHEAR_MODULUS[layer] = sm_real + sm_img * 1j
                LAME_CONSTANT[layer] = BULK_MODULUS[layer] - 2*SHEAR_MODULUS[layer]/3
    elif model=='Andrade':
        for layer, radius in enumerate(RADIUS):
            if radius <= Radius[0]:
                DENSITY[layer] = Rho[0]
                GRAVITY[layer] = 4/3*NP.pi*G_Constant * Rho[0] * radius
                mu_tmp = Mu[0]
                eta_tmp = Eta[0]
                BULK_MODULUS[layer] = Kappa[0]
                SHEAR_MODULUS[layer] = 0 + 0j
                LAME_CONSTANT[layer] = BULK_MODULUS[layer] - 2*SHEAR_MODULUS[layer]/3
            else:
                layer_material = LA.Find_MaterialLayer(radius.real, rr)
                DENSITY[layer] = Rho[layer_material]
                GRAVITY[layer] = 4/3*NP.pi*G_Constant*DENSITY[layer]*radius +  GRAVITY[layer-1]*(RADIUS[layer-1]/radius)**2 - 4/3*NP.pi*G_Constant*DENSITY[layer]*RADIUS[layer-1]**3/radius**2 
                mu_tmp = Mu[layer_material]
                eta_tmp = Eta[layer_material]
                BULK_MODULUS[layer] = Kappa[layer_material]
                Beta = mu_tmp**(Alpha-1) / eta_tmp**(Alpha) ##3.2e-13
                andrade_cos = Beta * SSP.gamma(1+Alpha) * Angular_Velocity**(-Alpha) * NP.cos(Alpha*NP.pi/2)
                andrade_sin = Beta * SSP.gamma(1+Alpha) * Angular_Velocity**(-Alpha) * NP.sin(Alpha*NP.pi/2)
                sm_real = mu_tmp * (1+mu_tmp*andrade_cos) / ((1+mu_tmp*andrade_cos)**2 + mu_tmp**2*(andrade_sin+1/(eta_tmp*Angular_Velocity))**2)
                sm_img = mu_tmp**2 * (andrade_sin+1/(eta_tmp*Angular_Velocity)) / ((1+mu_tmp*andrade_cos)**2 + mu_tmp**2*(andrade_sin+1/(eta_tmp*Angular_Velocity))**2)
                SHEAR_MODULUS[layer] = sm_real + sm_img * 1j
                LAME_CONSTANT[layer] = BULK_MODULUS[layer] - 2*SHEAR_MODULUS[layer]/3
    elif model=='none': ##zwb 20240830
        return 0.0, 0.0
    #----------------------------------------------------------------------------#
    #Generating Propagator Matrix based on Sabadini, Vermeersen & Cambiotti, 2016#
    #----------------------------------------------------------------------------#
    pgYk = NP.zeros((Layer_Number,6,6), dtype=complex) ## The propagator matrix (type: ndarray), product of k=0 to N-1, do not calculate Layer 0. 
    r = SY.Symbol('r', complex=True)                  ## The shape of ndarray can be >= 3 while the shape of matrix must be 2. 
    g = SY.Symbol('g', complex=True)                  ## Every element of pgYk should be converted to matrix before involved  
    mu = SY.Symbol('mu', complex=True) ## complex number
    rho = SY.Symbol('rho', complex=True)
    Y = SYM.Matrix([[r**3 / 7, r, 0, 1 / (2*r**2), 1 / r**4, 0],
                [5 * r**3 / 42, r / 2, 0, 0, -1 / (3*r**4), 0],
                [(-mu + g*rho*r) *r**2 / 7, 2*mu + g*rho*r, r**2 * rho, (-6*mu + g*rho*r) / (2*r**3), (-8*mu + g*rho*r) / r**5, rho / r**3],
                [8*mu * r**2 / 21, mu, 0, mu / (2*r**3), 8*mu / (3*r**5), 0],
                [0, 0, r**2, 0, 0, 1 / r**3],
                [4*SY.pi*G_Constant*rho*r**3 / 7, 4*SY.pi*G_Constant*rho*r, 5*r, 2*SY.pi*G_Constant*rho / r**2, 4*SY.pi*G_Constant*rho / r**4, 0]])
    D = SYM.diag(3 / r**3, 1 / r, -1 / r, 2*r**2, 3*r**4 / 7, r**3) / 5
    Ybar = SYM.Matrix([[rho * g * r / mu - 8, 16, -r / mu, 2 * r / mu, rho * r / mu, 0],
                    [-rho * g * r / mu + 6, -6, r / mu, 0, -rho * r / mu, 0],
                    [4 * SY.pi * G_Constant * rho, 0, 0, 0, 0, -1],
                    [rho * g * r / mu + 2, 6, -r / mu, -3 * r / mu, rho * r / mu, 0],
                    [-rho * g * r / mu + 1, -16, r / mu, 5 * r / mu, -rho * r / mu, 0],
                    [4 * SY.pi * G_Constant * rho * r, 0, 0, 0, 5, -r]])
    invY = D * Ybar ## Y.inv() works, too
    func_Yk = SY.lambdify((r, g, mu, rho), Y, modules='numpy') ## Convert Sympy matrices into Numpy ndarrays (i.e., matrices)
    func_invYk = SY.lambdify((r, g, mu, rho), invY, modules='numpy') ## Ditto
    for k, ite_tmp in enumerate(RADIUS[1:], start=1): ##enumerate will return a tuple, which is immutable. So we directly change the pgYk via index
        if k == 1:
            pgYk[k] = func_invYk(complex(RADIUS[0]), complex(GRAVITY[1]), complex(SHEAR_MODULUS[1]), complex(DENSITY[1]))
            continue
        else:
            Yk_tmp = func_Yk(complex(RADIUS[k-1]), complex(GRAVITY[k-1]), complex(SHEAR_MODULUS[k-1]), complex(DENSITY[k-1]))
            invYk_tmp = func_invYk(complex(RADIUS[k-1]), complex(GRAVITY[k]), complex(SHEAR_MODULUS[k]), complex(DENSITY[k]))
            pgYk[k] = invYk_tmp @ Yk_tmp @ pgYk[k-1] ## @ operator is matrix multiplication

    #--------------------------------------#
    #Calculating Interface Matrix#
    #CMB boundary condition based on#
    #Sabadini, Vermeersen & Cambiotti, 2016#
    #--------------------------------------#
    Ic = NP.array([[-Radius[0]**2 / complex(GRAVITY[0]), 0, 1], 
                [0, 1, 0],
                [0, 0, complex(GRAVITY[0]) * complex(DENSITY[0])],
                [0, 0, 0],
                [Radius[0]**2, 0, 0],
                [2 * Radius[0], 0, 4 * NP.pi* G_Constant * complex(DENSITY[0])]], dtype=complex) 

    #--------------------------------------#
    #Solving Cc, the constant vector#
    #Surface boundary condition based on#
    #Sabadini, Vermeersen & Cambiotti, 2016#
    #--------------------------------------#
    P1 = NP.diag(NP.array([0, 0, 1, 1, 0, 1], dtype=complex)) ##Selection Matrix
    Bs = NP.array([[0], [0], [5 / Radius[-1]]], dtype=complex) ##y3=0; y4=0; y6=5/Rs. ##test
    pgY = P1 @ func_Yk(complex(RADIUS[Layer_Number-1]), complex(GRAVITY[Layer_Number-1]), complex(SHEAR_MODULUS[Layer_Number-1]), complex(DENSITY[Layer_Number-1])) @ pgYk[Layer_Number-1] @ Ic
    Cc = NPL.solve(pgY[[2,3,5],:], Bs)

    #-----------------------------------------------------------------------------------------------------#
    #Calculating the vector y#
    #Note: the signs of y1-y4 obtained here are opposite to those in Tobie et al. (2005)#
    #This is because we define the outward radial direction as positive, i.e., tensile stress is positive#
    #while Tobie et al. (2005) consider the inward radia direction as positive, i.e., compress is positive#
    #NNote: In Tobie et al. (2005), y2 represented spheroidal radial stress#
    #and y3 represented spheroidal tangential displacement.#
    #-----------------------------------------------------------------------------------------------------#
    y1 = NP.zeros((Layer_Number,1), dtype=complex) ## spheroidal radial displacement
    y2 = NP.zeros((Layer_Number,1), dtype=complex) ## spheroidal tangential displacement
    y3 = NP.zeros((Layer_Number,1), dtype=complex) ## spheroidal radial stress
    y4 = NP.zeros((Layer_Number,1), dtype=complex) ## spheroidal tangential stress 
    y5 = NP.zeros((Layer_Number,1), dtype=complex) ## potential
    y6 = NP.zeros((Layer_Number,1), dtype=complex) ## potential stress
    for k, radius in enumerate(RADIUS):
        Yk = func_Yk(complex(RADIUS[k]), complex(GRAVITY[k]), complex(SHEAR_MODULUS[k]), complex(DENSITY[k]))
        yk = Yk @ pgYk[k] @ Ic @ Cc
        y1[k] = yk[0]
        y2[k] = yk[1]
        y3[k] = yk[2]
        y4[k] = yk[3]
        y5[k] = yk[4]
        y6[k] = yk[5]

    #---------------------------------------------------------#
    #Calculate H_mu and tidal dissipation rate per unit volume#
    #This method works for computating the radial distribution#
    #of the dissipation rate within any planetary interior#
    #according to Kervazo et al., 2021; Tobie et al., 2005#
    #---------------------------------------------------------#
    y_1 = SY.Symbol('y_1', complex=True)
    y_2 = SY.Symbol('y_2', complex=True)
    y_3 = SY.Symbol('y_3', complex=True)
    y_4 = SY.Symbol('y_4', complex=True)
    kappa = SY.Symbol('kappa', complex=True)
    dy_1 = (y_2 - (kappa-2*mu/3)*(2*y_1-6*y_3)/r) / (kappa + 4*mu/3) 
    H_mu_part1 = 4/3 * (r/SY.Abs(kappa+4*mu/3))**2 * (SY.Abs(y_2 - (kappa-2*mu/3)*(2*y_1 - 6*y_3)/r))**2
    H_mu_part2 = -4/3 * r * (dy_1.conjugate()*(2*y_1 - 6*y_3)).as_real_imag()[0]
    H_mu_part3 = 1/3 * (SY.Abs(2*y_1 - 6*y_3))**2 
    H_mu_part4 = 6 * r**2 * (SY.Abs(y_4))**2 / (SY.Abs(mu))**2 
    H_mu_part5 = 24 * (SY.Abs(y_3))**2 
    H_mu = H_mu_part1 + H_mu_part2 + H_mu_part3 + H_mu_part4 + H_mu_part5
    func_Hmu = SY.lambdify((y_1, y_2, y_3, y_4, mu, kappa, r), H_mu, modules='numpy')

    #------------------------------------------------------------#
    #Calculating the tidal dissipation rate#
    #according to Eq. (37) in Tobie et al. (2005)#
    #Note: our y1-y4 are opposite to those of Tobie et al. (2005)#
    #and y2-y3 should exchange #
    #------------------------------------------------------------#
    Hmu = NP.zeros((Layer_Number,1)) 
    h_tide = NP.zeros((Layer_Number,1)) ## tidal dissipation rate, W/m^3
    Volume_MIC = 4/3 * NP.pi * (rr[1]**3 - rr[0]**3)
    Volume_mantle = 4/3 * NP.pi * (rr[2]**3 - rr[1]**3)
    Power_tide_MIC = 0
    Power_tide_mantle = 0
    for k, radius in enumerate(RADIUS):
        if radius.real <= rr[0]:
            continue
        elif radius.real <= rr[1]:
            Hmuk = func_Hmu(complex(-y1[k]), complex(-y3[k]), complex(-y2[k]), complex(-y4[k]), 
                            complex(SHEAR_MODULUS[k]), complex(BULK_MODULUS[k]), complex(radius))
            Hmu[k] = Hmuk.real
            h_tide[k] = 21/10 * Angular_Velocity**5 * Rs**4 * Eccentricity**2 / radius**2 * Hmu[k] * SHEAR_MODULUS[k].imag
            Power_tide_MIC = Power_tide_MIC + 4/3*NP.pi * (radius**3 - RADIUS[k-1]**3) * h_tide[k]
        elif radius.real <= rr[2]:
            Hmuk = func_Hmu(complex(-y1[k]), complex(-y3[k]), complex(-y2[k]), complex(-y4[k]), 
                            complex(SHEAR_MODULUS[k]), complex(BULK_MODULUS[k]), complex(radius))
            Hmu[k] = Hmuk.real
            h_tide[k] = 21/10 * Angular_Velocity**5 * Rs**4 * Eccentricity**2 / radius**2 * Hmu[k] * SHEAR_MODULUS[k].imag
            Power_tide_mantle = Power_tide_mantle + 4/3*NP.pi * (radius**3 - RADIUS[k-1]**3) * h_tide[k]
        else:
            Hmuk = func_Hmu(complex(-y1[k]), complex(-y3[k]), complex(-y2[k]), complex(-y4[k]), 
                            complex(SHEAR_MODULUS[k]), complex(BULK_MODULUS[k]), complex(radius))
            Hmu[k] = Hmuk.real
            h_tide[k] = 21/10 * Angular_Velocity**5 * Rs**4 * Eccentricity**2 / radius**2 * Hmu[k] * SHEAR_MODULUS[k].imag

    return (Power_tide_MIC[0].real/Volume_MIC), (Power_tide_mantle[0].real/Volume_mantle)

def func_interpolate_Earth_Moon_distance(distance_model):
    #-----------------------------#
    #Basic Parameters for the Moon#
    #-----------------------------#
    Radius_Earth = 6371e3
    Distance_EM = 384400e3
    Time_unit = 365 * 24 * 60 * 60 * 1e9
    Distance_unit = Distance_EM / Radius_Earth
    if distance_model=='Webb1982':
        #----------------------------------------------#
        #calculate Webb, 1982#
        #https://github.com/trichter/archean_moon_orbit#
        #----------------------------------------------#
        models = NP.load('Data_From_Eulenfeld&Heubeck/data/orbit_models.npz')
        data_Webb = models['Webb 1982 curve d'] ## Time (Ga); Distance (normalized by present-day value)
        time_Webb = data_Webb[0]
        a_RE_Webb = data_Webb[1] * Distance_unit
        func_interp_a_RE = CSINTP(time_Webb, a_RE_Webb)

    elif distance_model=='Daher2021_PD':
        #---------------------------------------------------------------#
        #calculate Daher et al., 2021        #
        #data interpolation & derivation#
        #https://deepblue.lib.umich.edu/data/concern/data_sets/sj1392193#
        #---------------------------------------------------------------#
        #ocean basin geometries of present-day#
        data_PD = scipy.io.loadmat('Data_Daher2021/integration_results_use_Schindelegger_PD_experiments_ode45_nodeLF_v9.mat')
        time_PD = data_PD['timestep_vector_PD_ode45'].flatten() / (86400 * 365.25 * 10**9)
        a_RE_PD = data_PD['StateVector_PD_ode45'][:, 2] / 6378136
        func_interp_a_RE = CSINTP(time_PD, a_RE_PD)

    elif distance_model=='Tyler2021_40':
        #----------------------------------------------#
        #calculate Tyler et al., 2021#
        #https://github.com/trichter/archean_moon_orbit#
        #----------------------------------------------#
        models = NP.load('Data_From_Eulenfeld&Heubeck/data/orbit_models.npz')
        data_Tyler_T40 = models['Tyler 2021 T=40']
        #Best-fit values: h=2.3 km and T=40#
        a_RE_Tyler2021_40 =  data_Tyler_T40[1]*Distance_unit
        func_interp_a_RE = CSINTP(data_Tyler_T40[0], a_RE_Tyler2021_40)

    elif distance_model=='Farhat2022_mid':
        #------------------------------------#
        #calculate Farhat et al., 2022#
        #https://www.astrogeo.eu/?page_id=553#
        #------------------------------------#
        path_data_Farhat = ('Farhat_2022.csv')
        data_Farhat = NP.loadtxt(path_data_Farhat, delimiter=',', skiprows=1) 
        time_Farhat = data_Farhat[:,0] #Column 0: Ga; 
        a_RE_middle = data_Farhat[:,1] #Column 1: a(RE);
        func_interp_a_RE = CSINTP(time_Farhat, a_RE_middle)

    elif distance_model=='Farhat2022_min':
        #------------------------------------#
        #calculate Farhat et al., 2022#
        #https://www.astrogeo.eu/?page_id=553#
        #------------------------------------#
        path_data_Farhat = ('Farhat_2022.csv')
        data_Farhat = NP.loadtxt(path_data_Farhat, delimiter=',', skiprows=1) 
        time_Farhat = data_Farhat[:,0] #Column 0: Ga; 
        a_RE_min = data_Farhat[:,3] #Column 3: amin
        func_interp_a_RE = CSINTP(time_Farhat, a_RE_min)

    elif distance_model=='Farhat2022_max':
        #------------------------------------#
        #calculate Farhat et al., 2022#
        #https://www.astrogeo.eu/?page_id=553#
        #------------------------------------#
        path_data_Farhat = ('Farhat_2022.csv')
        data_Farhat = NP.loadtxt(path_data_Farhat, delimiter=',', skiprows=1) 
        time_Farhat = data_Farhat[:,0] #Column 0: Ga; 
        a_RE_max = data_Farhat[:,2] #Column 2: amax;
        func_interp_a_RE = CSINTP(time_Farhat, a_RE_max)

    return func_interp_a_RE