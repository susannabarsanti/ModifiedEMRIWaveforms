import numpy as np

from few.trajectory.ode import KerrEccEqFlux
from few.utils.geodesic import  get_separatrix

from multispline.spline import BicubicSpline

from pathlib import Path

# =====================================================================
# Scalar flux for circular equatorial Kerr orbits
# =====================================================================

class KerrCircEqFluxScalar(KerrEccEqFlux):
    @staticmethod
    # -----------------------------------------------------------------
    # Load scalar flux table from file and apply PN rescaling
    # -----------------------------------------------------------------
    def load_scalar_flux_table(filename):
        """
        Load a scalar flux table from a .dat file and apply a post-Newtonian rescaling.
        
        Parameters
        ----------
        filename : str or Path
            Path to the scalar flux table file.
            
        Returns
        -------
        z_col, u_col : arrays
            Grid coordinates used for interpolation.
        F_total_with_PN_rescaling : array
            Scalar flux (horizon + infinity) with PN rescaling applied.
        """
        
        # Load table data
        data = np.loadtxt(filename)
        a_col = data[:, 0]      # spin of the black hole
        p_col = data[:, 1]      # semi-latus rectum
        F_hor = data[:, 4]      # scalar flux at the horizon
        F_inf = data[:, 5]      # scalar flux at infinity
        z_col = data[:, 7]      # rescaled spin coordinate for interpolation
        u_col = data[:, 8]      # rescaled p coordinate for interpolation

        # Total scalar flux
        F_total_col = (F_hor + F_inf)

        # Leading-order post-Newtonian flux
        Epn = 1.0 / (12. * p_col**(4.))

        # Rescale the flux for interpolation
        F_total_with_PN_rescaling = (F_total_col - Epn) * p_col**(6.)
    
        return z_col, u_col, F_total_with_PN_rescaling


    # -----------------------------------------------------------------
    # Initialize object: load flux tables and precompute constants
    # -----------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        """
        Initialize the flux object. Precompute constants, separatrix, and table spin.
        """

        # Call the parent class constructor (KerrEccEqFlux)
        super().__init__(*args, **kwargs)


        # Tolerance for interpolation bounds
        self.coord_tol = 1e-12  # used for both z and u checks

        # Paths to flux table files
        BASE_DIR = Path(__file__).resolve().parent
        DATA_DIR = BASE_DIR / "data"
        data_file_A = DATA_DIR / "tabA_ScalarKerr.dat"
        data_file_B = DATA_DIR / "tabB_ScalarKerr.dat"
        
        # Load tables
        z_A, u_A, F_A = self.load_scalar_flux_table(data_file_A) 
        z_B, u_B, F_B = self.load_scalar_flux_table(data_file_B)

        # Get unique grid points for bicubic interpolation
        z_A_unique_vals = np.unique(z_A)    
        u_A_unique_vals = np.unique(u_A)   
        z_B_unique_vals = np.unique(z_B)    
        u_B_unique_vals = np.unique(u_B)   
        
        # Reshape flux arrays into 2D grids (z x u)
        F_A_grid = F_A.reshape((len(z_A_unique_vals), len(u_A_unique_vals)))
        F_B_grid = F_B.reshape((len(z_B_unique_vals), len(u_B_unique_vals)))
        
        # Build bicubic interpolators for the two regions of the flux table
        self.Fphi_interp_A = BicubicSpline(z_A_unique_vals, u_A_unique_vals, F_A_grid)
        self.Fphi_interp_B = BicubicSpline(z_B_unique_vals, u_B_unique_vals, F_B_grid)   
        
        # Precompute constants for interpolation
        self.amin = -999/1000
        self.amax = 999/1000
        self.chi_min = (1 - self.amax)**(1/3)
        self.chi_max = (1 - self.amin)**(1/3)

        self.delta_pAmin = 1e-3
        self.delta_pAmax = 9.0 + self.delta_pAmin
        self.Cp = self.delta_pAmax - 2 * self.delta_pAmin
        self.CDelta = np.log(self.delta_pAmax - self.delta_pAmin)
        self.delta_pBmin = self.delta_pAmax
        self.delta_pBmax = 200.0

    # -----------------------------------------------------------------
    # Compute scalar energy flux
    # -----------------------------------------------------------------
    def compute_Edot_phi(self, prograde, y):
        """
        Compute scalar energy flux Edot_phi for a circular orbit at radius p.
        """

        p = y[0]
        e = y[1]
        x = y[2]

        if (prograde == True):
            a_table = np.abs(self.a)

        else:
            a_table = -np.abs(self.a)



        if not np.isclose(e, 0.0):
            raise ValueError("Eccentricity e must be zero for circular orbits.")

        # Compute separatrix for table spin
        p_sep = get_separatrix(a_table, 0.0, 1.0) #a_table is positve if orbit is prograde and negative is orbit is retrograde, so to compute the separatrix we can just take x = 1.0. 

        
        # Precompute max p for interpolation regions
        pmaxA = self.delta_pAmax + p_sep
        pmaxB = self.delta_pBmax + p_sep

        if p <= p_sep:
            raise ValueError(f"Orbital radius p={p} must be > separatrix p_sep={p_sep}.")

        # Define normalized spin coordinate for interpolation
        z = ((1 - a_table)**(1/3) - self.chi_min) / (self.chi_max - self.chi_min)
        
        # Check if z is within valid interpolation range
        if z < - self.coord_tol or z > 1.0 + self.coord_tol:
            raise ValueError(f"Normalized spin coordinate z = {z} is outside valid range [{- self.coord_tol}, {1.0 + self.coord_tol}].")

        # Clip z to [0, 1] to ensure it stays within interpolation bounds
        z = np.clip(z, 0.0, 1.0)

        # Compute normalized orbital radius coordinate for interpolation
        #  Interpolate using table A
        if p < pmaxA:
            u = (np.log(p - p_sep + self.Cp) - self.CDelta) / np.log(2)
            
            if u < - self.coord_tol or u > 1.0 + self.coord_tol:
                raise ValueError(f"Normalized spin coordinate u = {u} is outside valid range [{- self.coord_tol}, {1.0 + self.coord_tol}].")
            
            # Clip u to [0, 1] to ensure it stays within interpolation bounds
            u = np.clip(u, 0.0, 1.0)
            
            F = self.Fphi_interp_A(z, u)   
        #  Interpolate using table B              
        else:
            u = (self.delta_pBmin**(-0.5) - (p - p_sep)**(-0.5)) / (self.delta_pBmin**(-0.5) - (pmaxB - p_sep)**(-0.5))

            if u < -self.coord_tol or u > 1.0 + self.coord_tol:
                raise ValueError(f"Normalized spin coordinate u = {u} is outside valid range [{-self.coord_tol}, {1.0 + self.coord_tol}].")
            
            # Clip u to [0, 1] to ensure it stays within interpolation bounds
            u = np.clip(u, 0.0, 1.0)

            F = self.Fphi_interp_B(z, u)           
        
        # Add leading-order post-Newtonian contribution 
        Epn = 1.0 / (12. * p**(4.))
        Edot_phi = (F * p**(- 6.) + Epn)

        return Edot_phi

    
    
    # -----------------------------------------------------------------
    # Modify ODE RHS with scalar flux
    # -----------------------------------------------------------------
    def modify_rhs(self, ydot, y):
        """
        Add scalar flux contribution to the RHS of orbital evolution ODE.
        """
        x = y[2]
    
        if (self.a * x >= 0.0):
            prograde = True
            a_evolution = np.abs(self.a)
        else:
            prograde = False
            a_evolution = -np.abs(self.a)

        # Compute scalar flux at current radius
        Edot_phi = self.compute_Edot_phi(prograde, y)

        # Scale by scalar charge squared if provided
        q_s2 = self.additional_args[0] if hasattr(self, "additional_args") and len(self.additional_args) > 0 else 0.0
        Edot_phi *= q_s2
        
        # Compute dE/dp using Kerr circular orbit formula
        sqrt_r = np.sqrt(y[0])
        r_32 = y[0]**1.5
        r2 = y[0]**2
        numerator = -3 * a_evolution**2 + 8 * a_evolution * sqrt_r + (-6 + y[0]) * y[0]
        term1 = 2 * a_evolution * y[0] + (-3 + y[0]) * r_32
        term2 = 2 * a_evolution * r_32 + (-3 + y[0]) * r2
        denominator = 2 * term1 * np.sqrt(term2)
        dE_dp = numerator/denominator
        
        # Update RHS of orbital evolution ODE with scalar flux contribution
        ydot[0] -= Edot_phi / dE_dp
