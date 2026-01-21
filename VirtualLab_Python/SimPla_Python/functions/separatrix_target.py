import numpy as np
from matplotlib.path import Path

class separatrix_target:
    def __init__(self):
        self.R_sep_target = None
        self.Z_sep_target = None
        self.inside = None
        self.touching_wall = 0

    def build_separatrix(self, separatrix_config, geo):
    
        if separatrix_config.method == 1:
            self.build_separatrix_m1(separatrix_config, geo)
            self.separatrix_outside_wall_v1(geo)
        # Add future methods here

        self.inside_separatrix(geo)
        return self

    def inside_separatrix(self, geo):
        Rg = geo.grid.Rg
        Zg = geo.grid.Zg
        points = np.vstack((Rg.ravel(), Zg.ravel())).T

        polygon = np.column_stack((self.R_sep_target, self.Z_sep_target))
        path = Path(polygon)
        inside = path.contains_points(points).reshape(Rg.shape)

        self.inside = inside

    def sep_operators(self, geo):
        Rg = geo.grid.Rg
        Zg = geo.grid.Zg
        R_sep = self.R_sep_target
        Z_sep = self.Z_sep_target

        ind_sep = []
        for r_s, z_s in zip(R_sep, Z_sep):
            dist2 = (Rg - r_s)**2 + (Zg - z_s)**2
            ind = np.argmin(dist2)
            ind_sep.append(ind)

        ind_sep = np.unique(ind_sep)
        M_sep = np.zeros((len(ind_sep), Rg.size))
        V_sep = np.zeros((len(ind_sep), 1))

        for i, idx in enumerate(ind_sep):
            M_sep[i, idx] = 1

        return M_sep, V_sep, ind_sep

    def build_separatrix_m1(self, config, geo):
        
        p = config
        
        k1, k2 = p.k1, p.k2
        d1, d2 = p.d1, p.d2
        gamma_n1, gamma_n2 = p.gamma_n_1, p.gamma_n_2
        gamma_p1, gamma_p2 = p.gamma_p_1, p.gamma_p_2
    
        R0 = p.R0
        Z0 = p.Z0
        a = p.a
        
        def inner_segment(k, d, gamma, upper=True):

            sign = 1 if upper else -1
            tn = (1 - d) / k * np.tan(gamma)
    
            if tn < 0.5:
                alpha0 = -(d - (1 + d) * tn) / (1 - 2 * tn)
                alpha = (1 - d) * (1 - tn) / (1 - 2 * tn)
                beta = k * (1 - tn) / np.sqrt(1 - 2 * tn)
                thetax = np.arcsin(np.sqrt(1 - 2 * tn) / (1 - tn))
                theta = np.linspace(0, thetax, 100) * sign
                R = R0 + a * (alpha0 - alpha * np.cos(theta))
                Z = a * beta * np.sin(theta)
            elif tn == 0.5:
                zeta = np.linspace(0, k, 100) * sign
                xi = -1 + (1 + d) / k**2 * zeta**2
                R = xi * a + R0
                Z = zeta * a
            elif tn > 0.5 and tn < np.inf:
                alpha0 = -((1 + d) * tn - d) / (2 * tn - 1)
                alpha = (1 - d) * (1 - tn) / (2 * tn - 1)
                beta = k * (1 - tn) / np.sqrt(2 * tn - 1)
                phix = np.arcsinh(np.sqrt(2 * tn - 1) / (1 - tn))
                phi = np.linspace(0, phix, 100) * sign
                R = R0 + a * (alpha0 + alpha * np.cosh(phi))
                Z = a * beta * np.sinh(phi)
            elif tn == 1:
                zeta = np.linspace(0, k, 100) * sign
                xi = -1 + (1 - d) / k * zeta
                R = xi * a + R0
                Z = zeta * a
            else:
                R, Z = np.array([]), np.array([])
                
            return R, Z
            
        def outer_segment(k, d, gamma, upper=True):
            sign = 1 if upper else -1
            tp = (1 + d) / k * np.tan(gamma)
            
            if tp < 0.5:
                alpha0 = -(d + (1 - d) * tp) / (1 - 2 * tp)
                alpha = (1 + d) * (1 - tp) / (1 - 2 * tp)
                beta = k * (1 - tp) / np.sqrt(1 - 2 * tp)
                thetax = np.arcsin(np.sqrt(1 - 2 * tp) / (1 - tp))
                theta = np.linspace(0, thetax, 100) * sign
                R = R0 + a * (alpha0 + alpha * np.cos(theta))
                Z = a * beta * np.sin(theta)
            elif tp == 0.5:
                zeta = np.linspace(0, k, 100) * sign
                xi = -1 - (1 + d) / k**2 * zeta**2
                R = xi * a + R0
                Z = zeta * a
            elif tp > 0.5 and tp < 1:
                alpha0 = ((1 - d) * tp + d) / (2 * tp - 1)
                alpha = -(1 + d) * (1 - tp) / (2 * tp - 1)
                beta = k * (1 - tp) / np.sqrt(2 * tp - 1)
                phix = np.arcsinh(np.sqrt(2 * tp - 1) / (1 - tp))
                phi = np.linspace(0, phix, 100) * sign
                R = R0 + a * (alpha0 + alpha * np.cosh(phi))
                Z = a * beta * np.sinh(phi)
            elif tp == 1:
                zeta = np.linspace(0, k, 100) * sign
                xi = 1 - (1 + d) / k * zeta
                R = xi * a + R0
                Z = zeta * a
            else:
                R, Z = np.array([]), np.array([])
            return R, Z
        
        # Calcolo sezioni
        Rnu, Znu = inner_segment(k1, d1, gamma_n1, upper=True)
        Rnl, Znl = inner_segment(k2, d2, gamma_n2, upper=False)
        Rpu, Zpu = outer_segment(k1, d1, gamma_p1, upper=True)
        Rpl, Zpl = outer_segment(k2, d2, gamma_p2, upper=False)
        
        # Composizione finale
        self.R_sep_target = np.concatenate([Rnu, Rpu[::-1], Rpl, Rnl[::-1]])
        self.Z_sep_target = np.concatenate([Znu, Zpu[::-1], Zpl, Znl[::-1]]) + Z0
            
    def separatrix_outside_wall_v1(self, geo):
        """
        Correct target separatrix if it crosses the wall
        """

        R_sep = np.asarray(self.R_sep_target)
        Z_sep = np.asarray(self.Z_sep_target)

        # Polygon wall
        wall_poly = np.column_stack((geo.wall.R, geo.wall.Z))
        sep_points = np.column_stack((R_sep, Z_sep))

        inside = Path(wall_poly).contains_points(sep_points)

        # Check if separatrix touches wall
        if np.sum(inside) < len(inside):
            self.touching_wall = 1
        else:
            self.touching_wall = 0

        # Keep only points inside the wall
        self.R_sep_target = R_sep[inside]
        self.Z_sep_target = Z_sep[inside]
           
