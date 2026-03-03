/* rigid_cavity_model.c
 *
 * BENG207_2: Rigid Cavity Resonator Model
 *
 * Lateral standing wave between rigid sidewalls for iPSC organoid formation.
 * Organoid spacing ~1 mm matches CI electrode spacing (~30 mm cochlea / 22 ch).
 * Operating frequency ~741 kHz (lambda/2 = 1 mm in water).
 *
 * Stack (bottom to top, vertical energy coupling):
 *   1. LiNbO3 transducer (source)
 *   2. Ultrasound gel couplant (Kelvin-Voigt)
 *   3. Borosilicate glass superstrate
 *   4. Fluid channel (water) — lateral standing wave forms HERE
 *
 * The standing wave propagates LATERALLY (along x, parallel to glass)
 * between two rigid sidewalls, NOT vertically through the stack.
 *
 * See README file and BENG207_model_2_schematic.svg for geometry.
 *
 * Part A: Lateral cavity — pressure field, contrast, radiation force
 * Part B: Vertical coupling — impedance transformer (from BENG207_1)
 *
 * Compile (Linux/Mac):
 *   gcc -std=c99 -O2 -o cavity_model rigid_cavity_model.c -lm
 *
 * Compile (Windows MSVC):
 *   cl /O2 rigid_cavity_model.c /link /out:cavity_model.exe
 *
 * Compile (Windows MinGW):
 *   gcc -O2 -o cavity_model.exe rigid_cavity_model.c -lm
 *
 * Run:
 *   ./cavity_model        (Linux/Mac)
 *   cavity_model.exe      (Windows)
 *   CSVs appear in PLOT_CSV/ subfolder next to the executable.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================
 * Platform-specific directory creation
 * ============================================================ */
#ifdef _WIN32
  #include <direct.h>
  #define MKDIR(path) _mkdir(path)
#else
  #include <sys/stat.h>
  #include <errno.h>
  #define MKDIR(path) mkdir(path, 0755)
#endif

/* ============================================================
 * Complex number helpers (MSVC lacks <complex.h>)
 * ============================================================ */

typedef struct { double re; double im; } Cpx;

static Cpx cpx(double re, double im)
{
    Cpx z; z.re = re; z.im = im; return z;
}

static Cpx cpx_add(Cpx a, Cpx b)
{
    return cpx(a.re + b.re, a.im + b.im);
}

static Cpx cpx_sub(Cpx a, Cpx b)
{
    return cpx(a.re - b.re, a.im - b.im);
}

static Cpx cpx_mul(Cpx a, Cpx b)
{
    return cpx(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

static Cpx cpx_div(Cpx a, Cpx b)
{
    double d = b.re*b.re + b.im*b.im;
    return cpx((a.re*b.re + a.im*b.im) / d,
               (a.im*b.re - a.re*b.im) / d);
}

static Cpx cpx_scale(double s, Cpx a)
{
    return cpx(s * a.re, s * a.im);
}

static double cpx_abs(Cpx a)
{
    return sqrt(a.re*a.re + a.im*a.im);
}

static Cpx cpx_exp(Cpx a)
{
    double e = exp(a.re);
    return cpx(e * cos(a.im), e * sin(a.im));
}

static Cpx cpx_conj(Cpx a)
{
    return cpx(a.re, -a.im);
}

static Cpx cpx_sin(Cpx a)
{
    return cpx(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}

static Cpx cpx_cos(Cpx a)
{
    return cpx(cos(a.re) * cosh(a.im), -sin(a.re) * sinh(a.im));
}

static Cpx cpx_sqrt(Cpx a)
{
    double r = sqrt(cpx_abs(a));
    double theta = atan2(a.im, a.re) / 2.0;
    return cpx(r * cos(theta), r * sin(theta));
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Output directory — relative to executable */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_Echo_Chamber_AMN/PLOT_CSV"

/* ============================================================
 * Material parameters
 * ============================================================ */
/* Parameters have been changed to reflect what we currently use */

/* Fluid channel (water) — where the lateral standing wave forms */
static const double rho_f = 1000.0;     /* kg/m^3 */
static const double c_f   = 1500;     /* m/s */

/* Fluid attenuation model: alpha = alpha_coeff * f^2
 * For water at 25C: alpha_coeff ~ 25e-15 Np/m/Hz^2
 * At 741 kHz: alpha ~ 25e-15 * (741e3)^2 ~ 1.4e-5 Np/m (negligible)
 * At 20 MHz:  alpha ~ 25e-15 * (20e6)^2  ~ 0.01 Np/m              */
static const double alpha_coeff = 25.0e-15;  /* Np/m/Hz^2 */

/* Wall materials */
/* Glass*/
static const double rho_glass = 2210.0;     /* kg/m^3 */
static const double c_glass   = 5700.0;     /* m/s */


static const double rho_si    = 2330.0;     /* kg/m^3 */
static const double c_si      = 8433.0;     /* m/s (longitudinal) */

/* ============================================================
 * BEGIN ADDITION: Gor'kov acoustic radiation force parameters
 *
 * The Gor'kov potential for a small sphere in a standing wave:
 *   U = (4/3)*pi*a^3 * [ f1/(2*rho_f*c_f^2) * <p^2>
 *                       - 3*f2*rho_f/(4)      * <v^2> ]
 *
 * where f1 = 1 - kappa_p/kappa_f   (monopole coefficient)
 *       f2 = 2*(rho_p - rho_f) / (2*rho_p + rho_f)  (dipole)
 *
 * For a 1-D standing wave p = p0*cos(kx), the radiation force
 * simplifies to:
 *   F_rad = 4*pi*a^3 * k * Phi * <E_ac> * sin(2kx)
 *
 * where Phi = f1/3 + f2/2  is the acoustic contrast factor,
 *       <E_ac> = p0^2 / (4*rho_f*c_f^2)  is the energy density.
 *
 * Phi > 0  =>  particles migrate to pressure NODES  (our case)
 * Phi < 0  =>  particles migrate to pressure ANTINODES
 *
 * References:
 *   Gor'kov, Sov. Phys. Doklady 6, 773 (1962)
 *   Bruus, Lab Chip 12, 1014 (2012)
 *   Cohen et al., Sci. Rep. 10, 4932 (2020) — neuronal patterning
 * ============================================================ */

/* Particle (iPSC organoid / neural aggregate) properties */
static const double rho_p = 1050.0;   /* kg/m^3  — cell density, typical
                                        * range 1040-1080 for mammalian cells.
                                        * Cohen et al. used DRG neurons
                                        * (~1050 kg/m^3). Adjust if measured. */

static const double c_p   = 1550.0;   /* m/s — longitudinal sound speed in
                                        * cell aggregates. Typical range
                                        * 1520-1600 m/s for soft tissue.
                                        * Higher than water due to protein
                                        * content. Adjust if measured. */

static const double a_organoid = 100.0e-6;  /* m — organoid radius (~200 um
                                              * diameter). This is the nominal
                                              * radius for force calculation;
                                              * actual organoids vary. */

/* Derived Gor'kov quantities (set in main) */
static double f1_monopole;  /* = 1 - kappa_p/kappa_f = 1 - (rho_f*c_f^2)/(rho_p*c_p^2) */
static double f2_dipole;    /* = 2*(rho_p - rho_f) / (2*rho_p + rho_f) */
static double Phi_contrast; /* = f1/3 + f2/2  — the acoustic contrast factor */

/* END ADDITION: Gor'kov parameters */

/* Vertical stack (for coupling sub-model, Part B) */
static const double rho_1    = 4647.0;      /* LiNbO3 kg/m^3 */
static const double c_1      = 6570.0;      /* LiNbO3 m/s */

static const double rho_2    = 1020.0;      /* couplant kg/m^3 */
static const double c_2      = 1500.0;      /* couplant m/s */

static const double rho_3    = 2230.0;      /* glass kg/m^3 */
static const double c_3      = 5640.0;      /* glass m/s */

/* Derived impedances (set in main) */
static double Z_f, Z_glass, Z_si, Z_1, Z_2, Z_3;

/* ============================================================
 * Utilities
 * ============================================================ */

static void linspace(double start, double end, int n, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = start + (end - start) * i / (n - 1);
}

static double dmax(double a, double b)
{
    return (a > b) ? a : b;
}

static void build_path(char *buf, int bufsize, const char *filename)
{
#ifdef _WIN32
    _snprintf(buf, bufsize, "%s\\%s", OUTPUT_DIR, filename);
#else
    snprintf(buf, bufsize, "%s/%s", OUTPUT_DIR, filename);
#endif
}

/* ============================================================
 * PART A: Lateral standing wave in rigid cavity
 *
 * p(x) = p0 * exp(-(alpha - ik)x) + r * p0 * exp(-(alpha - ik)(2L - x))
 *
 * where alpha = attenuation, k = 2*pi*f/c_f, r = reflection coeff,
 * L = channel length, x = position along channel.
 * ============================================================ */

/* Pressure at position x in the cavity */
static Cpx pressure_at(double x, double L, double freq,
                        double r_abs, double alpha)
{
    double k = 2.0 * M_PI * freq / c_f;
    Cpx gamma = cpx(-alpha, k);       /* propagation constant */
    Cpx r = cpx(r_abs, 0);            /* real reflection coeff */

    /* Forward wave: exp(-gamma * x) */
    Cpx fwd = cpx_exp(cpx_scale(-1.0, cpx_scale(x, gamma)));

    /* Reflected wave: r * exp(-gamma * (2L - x)) */
    Cpx ref = cpx_mul(r, cpx_exp(cpx_scale(-1.0, cpx_scale(2.0*L - x, gamma))));

    return cpx_add(fwd, ref);
}

/* Pressure squared (time-averaged): |p(x)|^2 */
static double pressure_sq(double x, double L, double freq,
                           double r_abs, double alpha)
{
    Cpx p = pressure_at(x, L, freq, r_abs, alpha);
    return p.re*p.re + p.im*p.im;
}

/* Acoustic radiation force proportional to -d/dx[<p^2>]
 * Computed by central finite difference */
static double radiation_force(double x, double L, double freq,
                               double r_abs, double alpha)
{
    double dx = L * 1e-5;  /* small step */
    double p2_plus, p2_minus;

    if (x - dx < 0) {
        p2_plus  = pressure_sq(x + dx, L, freq, r_abs, alpha);
        p2_minus = pressure_sq(x, L, freq, r_abs, alpha);
        return -(p2_plus - p2_minus) / dx;
    }
    if (x + dx > L) {
        p2_plus  = pressure_sq(x, L, freq, r_abs, alpha);
        p2_minus = pressure_sq(x - dx, L, freq, r_abs, alpha);
        return -(p2_plus - p2_minus) / dx;
    }

    p2_plus  = pressure_sq(x + dx, L, freq, r_abs, alpha);
    p2_minus = pressure_sq(x - dx, L, freq, r_abs, alpha);
    return -(p2_plus - p2_minus) / (2.0 * dx);
}

/* ============================================================
 * BEGIN ADDITION: Gor'kov radiation force (physical units)
 *
 * Converts the dimensionless pressure-gradient force into
 * absolute Newtons using the Gor'kov prefactor:
 *
 *   F_gorkov = Phi * (pi * a^3) / (3 * rho_f * c_f^2)
 *              * (-d/dx[<p^2>])
 *
 * This is derived from the Gor'kov potential gradient.
 * The factor (pi*a^3)/(3*rho_f*c_f^2) has units of
 * m^3 / (kg/m^3 * m^2/s^2) = m^3 * s^2 / kg = N / (Pa^2/m).
 *
 * To get force in Newtons, you still need to scale the
 * pressure field by the actual pressure amplitude p0 (Pa).
 * The output here assumes p0 = 1 Pa (i.e., the p field is
 * normalized). Multiply by p0^2 for actual force.
 *
 * Returns force in Newtons (for p0 = 1 Pa normalization).
 * Positive = toward nearest pressure node (trapping).
 * ============================================================ */
static double gorkov_force(double x, double L, double freq,
                            double r_abs, double alpha)
{
    /* Raw pressure-gradient force: -d/dx[<p^2>] in Pa^2/m */
    double F_raw = radiation_force(x, L, freq, r_abs, alpha);

    /* Gor'kov prefactor: Phi * pi * a^3 / (3 * rho_f * c_f^2) */
    double a3 = a_organoid * a_organoid * a_organoid;
    double prefactor = Phi_contrast * M_PI * a3
                       / (3.0 * rho_f * c_f * c_f);

    return prefactor * F_raw;
}
/* END ADDITION: Gor'kov radiation force */

/* Standing wave ratio at position x */
static double sw_ratio(double x, double L, double r_abs, double alpha)
{
    return r_abs * exp(-2.0 * alpha * (L - x));
}

/* Contrast from standing wave ratio */
static double contrast_from_swr(double S)
{
    if (S >= 0.9999) return 1e10;  /* avoid division by zero */
    return ((1.0 + S) / (1.0 - S)) * ((1.0 + S) / (1.0 - S));
}

/* Reflection coefficient from impedance mismatch */
static double reflection_coeff(double Z_wall, double Z_fluid)
{
    return fabs((Z_wall - Z_fluid) / (Z_wall + Z_fluid));
}

/* ============================================================
 * PART B: Vertical coupling (impedance transformer)
 *
 * Z_in(d) = Z2 * (Z3 + i*Z2*tan(k2*d)) / (Z2 + i*Z3*tan(k2*d))
 * T_power = 1 - |r_vert|^2
 * where r_vert = (Z_in - Z1) / (Z_in + Z1)
 * ============================================================ */

static double coupling_T_power(double freq, double d_couplant)
{
    double omega = 2.0 * M_PI * freq;
    double k2 = omega / c_2;
    double k2d = k2 * d_couplant;
    Cpx tan_k2d, Z_in, r_vert;
    Cpx numerator, denominator;
    double r_abs_sq;

    /* tan(k2d) — real since couplant is lossless in this sub-model */
    tan_k2d = cpx(tan(k2d), 0);

    /* Z_in = Z2 * (Z3 + i*Z2*tan(k2d)) / (Z2 + i*Z3*tan(k2d)) */
    numerator = cpx_add(cpx(Z_3, 0), cpx_mul(cpx(0, Z_2), tan_k2d));
    denominator = cpx_add(cpx(Z_2, 0), cpx_mul(cpx(0, Z_3), tan_k2d));
    Z_in = cpx_mul(cpx(Z_2, 0), cpx_div(numerator, denominator));

    /* r_vert = (Z_in - Z1) / (Z_in + Z1) */
    r_vert = cpx_div(cpx_sub(Z_in, cpx(Z_1, 0)),
                     cpx_add(Z_in, cpx(Z_1, 0)));

    r_abs_sq = r_vert.re*r_vert.re + r_vert.im*r_vert.im;
    return 1.0 - r_abs_sq;
}

/* Return Z_in as complex (for Fig 5 output) */
static Cpx coupling_Z_in(double freq, double d_couplant)
{
    double omega = 2.0 * M_PI * freq;
    double k2 = omega / c_2;
    double k2d = k2 * d_couplant;
    Cpx tan_k2d, numerator, denominator;

    tan_k2d = cpx(tan(k2d), 0);
    numerator = cpx_add(cpx(Z_3, 0), cpx_mul(cpx(0, Z_2), tan_k2d));
    denominator = cpx_add(cpx(Z_2, 0), cpx_mul(cpx(0, Z_3), tan_k2d));
    return cpx_mul(cpx(Z_2, 0), cpx_div(numerator, denominator));
}

/* ============================================================
 * OUTPUT 1: Lateral pressure field |p(x)|
 *
 * Three channel lengths (5, 10, 20 mm) at resonance,
 * for glass and silicon sidewalls.
 * ============================================================ */

static void output_pressure_field(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    double r_glass, r_silicon;
    int Nx = 2000;
    int iL, ix;

    build_path(fname, sizeof(fname), "fig1_pressure_field.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    fprintf(fp, "L_mm,x_mm,p_abs_glass,p_abs_silicon,p_sq_glass,p_sq_silicon\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        /* Choose mode number for ~1 mm spacing */
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;

        for (ix = 0; ix <= Nx; ix++) {
            double x = L * (double)ix / Nx;
            Cpx p_g = pressure_at(x, L, freq, r_glass, alpha);
            Cpx p_s = pressure_at(x, L, freq, r_silicon, alpha);

            fprintf(fp, "%.4f,%.6f,%.8f,%.8f,%.8f,%.8f\n",
                    L*1e3, x*1e3,
                    cpx_abs(p_g), cpx_abs(p_s),
                    p_g.re*p_g.re + p_g.im*p_g.im,
                    p_s.re*p_s.re + p_s.im*p_s.im);
        }
    }

    fclose(fp);
    printf("  Written: %s (3 lengths x %d positions)\n", fname, Nx+1);
}

/* ============================================================
 * OUTPUT 2: Standing wave contrast map
 *
 * 2D heatmap: contrast (dB) at channel midpoint
 * vs L (1-25 mm) and |r| (0.5-0.99)
 * ============================================================ */

static void output_contrast_map(void)
{
    char fname[256];
    FILE *fp;
    int NL = 200, Nr = 200;
    double L_arr[200], r_arr[200];
    int iL, ir;

    build_path(fname, sizeof(fname), "fig2_contrast_map.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e-3, 25e-3, NL, L_arr);
    linspace(0.50, 0.99, Nr, r_arr);

    fprintf(fp, "L_mm,r_abs,contrast_dB,contrast_linear,n_mode,S_mid\n");

    for (ir = 0; ir < Nr; ir++) {
        for (iL = 0; iL < NL; iL++) {
            double L = L_arr[iL];
            int n_mode = (int)(L / 1.0e-3 + 0.5);
            double freq = n_mode * c_f / (2.0 * L);
            double alpha = alpha_coeff * freq * freq;
            double S_mid = sw_ratio(L/2.0, L, r_arr[ir], alpha);
            double C = contrast_from_swr(S_mid);
            double C_dB = 10.0 * log10(dmax(C, 1e-15));

            fprintf(fp, "%.4f,%.4f,%.4f,%.4f,%d,%.6f\n",
                    L*1e3, r_arr[ir], C_dB, C, n_mode, S_mid);
        }
    }

    fclose(fp);
    printf("  Written: %s (%dx%d)\n", fname, NL, Nr);
}

/* ============================================================
 * OUTPUT 3: Resonant mode structure
 *
 * Mode frequencies f_n vs channel length,
 * number of modes in transducer bandwidth.
 * ============================================================ */

static void output_mode_structure(void)
{
    char fname[256];
    FILE *fp;
    int NL = 500;
    double L_arr[500];
    double bw_center = 741.5e3;  /* target frequency */
    double bw_half = 50.0e3;     /* assumed BW +/- 50 kHz */
    int iL;

    build_path(fname, sizeof(fname), "fig3_mode_structure.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e-3, 25e-3, NL, L_arr);

    fprintf(fp, "L_mm,n_target,f_target_kHz,delta_f_kHz,"
                "n_modes_in_bw,node_spacing_mm,organoid_count\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        double delta_f = c_f / (2.0 * L);  /* mode spacing */

        /* Mode number for ~1 mm node spacing */
        int n_target = (int)(L / 1.0e-3 + 0.5);
        if (n_target < 1) n_target = 1;

        double f_target = n_target * delta_f;
        double node_spacing = L / (double)n_target;

        /* Count modes within bandwidth */
        int n_lo = (int)ceil((bw_center - bw_half) / delta_f);
        int n_hi = (int)floor((bw_center + bw_half) / delta_f);
        int n_modes_in_bw = 0;
        if (n_hi >= n_lo && n_lo >= 1) n_modes_in_bw = n_hi - n_lo + 1;

        fprintf(fp, "%.4f,%d,%.4f,%.4f,%d,%.6f,%d\n",
                L*1e3, n_target, f_target/1e3, delta_f/1e3,
                n_modes_in_bw, node_spacing*1e3, n_target);
    }

    fclose(fp);
    printf("  Written: %s (%d lengths)\n", fname, NL);
}

/* ============================================================
 * OUTPUT 4: Nodal uniformity profile
 *
 * Radiation force magnitude at each node along the channel,
 * normalized to the strongest node.
 * Glass vs silicon, L = 5, 10, 20 mm.
 * ============================================================ */

static void output_nodal_uniformity(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    double r_glass, r_silicon;
    int iL, in;

    build_path(fname, sizeof(fname), "fig4_nodal_uniformity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    fprintf(fp, "L_mm,node_index,x_mm,"
                "p_sq_glass,p_sq_silicon,"
                "F_glass,F_silicon,"
                "F_glass_norm,F_silicon_norm,"
                /* BEGIN ADDITION: Gor'kov force columns (in Newtons, p0=1 Pa) */
                "F_gorkov_glass_N,F_gorkov_silicon_N\n");
                /* END ADDITION */

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;
        double node_spacing = L / (double)n_mode;
        double F_max_g = 0, F_max_s = 0;
        double *F_g, *F_s, *p2_g, *p2_s, *x_nodes;
        int i;

        F_g = (double *)malloc(n_mode * sizeof(double));
        F_s = (double *)malloc(n_mode * sizeof(double));
        p2_g = (double *)malloc(n_mode * sizeof(double));
        p2_s = (double *)malloc(n_mode * sizeof(double));
        x_nodes = (double *)malloc(n_mode * sizeof(double));

        /* Compute force at each node */
        for (in = 0; in < n_mode; in++) {
            double x = (in + 0.5) * node_spacing;  /* node centers */
            x_nodes[in] = x;
            p2_g[in] = pressure_sq(x, L, freq, r_glass, alpha);
            p2_s[in] = pressure_sq(x, L, freq, r_silicon, alpha);

            /* Force magnitude (absolute value of gradient) */
            F_g[in] = fabs(radiation_force(x, L, freq, r_glass, alpha));
            F_s[in] = fabs(radiation_force(x, L, freq, r_silicon, alpha));

            if (F_g[in] > F_max_g) F_max_g = F_g[in];
            if (F_s[in] > F_max_s) F_max_s = F_s[in];
        }

        /* Write normalized */
        for (i = 0; i < n_mode; i++) {
            /* BEGIN ADDITION: compute Gor'kov force at each node */
            double Fg_gorkov = fabs(gorkov_force(x_nodes[i], L, freq,
                                                  r_glass, alpha));
            double Fs_gorkov = fabs(gorkov_force(x_nodes[i], L, freq,
                                                  r_silicon, alpha));
            /* END ADDITION */

            fprintf(fp, "%.4f,%d,%.6f,%.8f,%.8f,%.8e,%.8e,%.6f,%.6f,"
                        /* BEGIN ADDITION: Gor'kov columns */
                        "%.8e,%.8e\n",
                        /* END ADDITION */
                    L*1e3, i, x_nodes[i]*1e3,
                    p2_g[i], p2_s[i],
                    F_g[i], F_s[i],
                    (F_max_g > 0) ? F_g[i]/F_max_g : 0,
                    (F_max_s > 0) ? F_s[i]/F_max_s : 0,
                    /* BEGIN ADDITION: Gor'kov values */
                    Fg_gorkov, Fs_gorkov);
                    /* END ADDITION */
        }

        free(F_g); free(F_s); free(p2_g); free(p2_s); free(x_nodes);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 5: Coupling layer characterization
 *
 * Z_in (real, imag, mag) and T_power vs couplant thickness
 * at 741 kHz.
 * ============================================================ */

static void output_coupling_layer(void)
{
    char fname[256];
    FILE *fp;
    int Nd = 2000;
    double d_arr[2000];
    double freq = 741.5e3;
    double lambda2 = c_2 / freq;  /* ~2.02 mm at 741 kHz */
    int i;

    build_path(fname, sizeof(fname), "fig5_coupling_layer.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(0.1e-6, 200e-6, Nd, d_arr);

    fprintf(fp, "d_um,d_over_lambda,Zin_re_MRayl,Zin_im_MRayl,"
                "Zin_mag_MRayl,T_power,T_power_dB\n");

    for (i = 0; i < Nd; i++) {
        Cpx Zin = coupling_Z_in(freq, d_arr[i]);
        double T = coupling_T_power(freq, d_arr[i]);
        double T_dB = 10.0 * log10(dmax(T, 1e-15));

        fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.6f,%.8f,%.4f\n",
                d_arr[i]*1e6, d_arr[i]/lambda2,
                Zin.re/1e6, Zin.im/1e6, cpx_abs(Zin)/1e6,
                T, T_dB);
    }

    fclose(fp);
    printf("  Written: %s (%d thicknesses, lambda2=%.3f mm at %.0f kHz)\n",
           fname, Nd, lambda2*1e3, freq/1e3);
}

/* ============================================================
 * OUTPUT 6: Frequency sensitivity
 *
 * Organoid spacing error and contrast loss vs frequency
 * detuning for different channel lengths.
 * ============================================================ */

static void output_frequency_sensitivity(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 15e-3, 20e-3};
    int NL = 4;
    int Ndf = 500;
    double df_arr[500];  /* fractional detuning */
    double r_wall;
    int iL, idf;

    build_path(fname, sizeof(fname), "fig6_freq_sensitivity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_wall = reflection_coeff(Z_glass, Z_f);
    linspace(-0.02, 0.02, Ndf, df_arr);  /* +/- 2% detuning */

    fprintf(fp, "L_mm,df_frac,df_kHz,f_actual_kHz,"
                "spacing_error_um,contrast_mid_dB,p_sq_mid\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double f_res = n_mode * c_f / (2.0 * L);

        for (idf = 0; idf < Ndf; idf++) {
            double f_actual = f_res * (1.0 + df_arr[idf]);
            double alpha = alpha_coeff * f_actual * f_actual;
            double lambda_actual = c_f / f_actual;
            double actual_spacing = lambda_actual / 2.0;
            double spacing_error = (actual_spacing - 1.0e-3) * 1e6;  /* um */

            /* Contrast at midpoint */
            double S_mid = sw_ratio(L/2.0, L, r_wall, alpha);
            double C_mid = contrast_from_swr(S_mid);
            double C_dB = 10.0 * log10(dmax(C_mid, 1e-15));

            /* Pressure squared at midpoint */
            double p2_mid = pressure_sq(L/2.0, L, f_actual, r_wall, alpha);

            fprintf(fp, "%.4f,%.6f,%.4f,%.4f,%.4f,%.4f,%.8f\n",
                    L*1e3, df_arr[idf], (f_actual - f_res)/1e3,
                    f_actual/1e3, spacing_error, C_dB, p2_mid);
        }
    }

    fclose(fp);
    printf("  Written: %s (%d lengths x %d detunings)\n", fname, NL, Ndf);
}

/* ============================================================
 * OUTPUT 7: Combined design space
 *
 * Organoid count, minimum nodal force, Q-factor vs channel
 * length for glass and silicon sidewalls.
 * ============================================================ */

static void output_design_space(void)
{
    char fname[256];
    FILE *fp;
    int NL = 500;
    double L_arr[500];
    double r_glass, r_silicon;
    int iL;

    build_path(fname, sizeof(fname), "fig7_design_space.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);
    linspace(1e-3, 25e-3, NL, L_arr);

    fprintf(fp, "L_mm,n_organoids,freq_kHz,"
                "Q_glass,Q_silicon,"
                "contrast_mid_dB_glass,contrast_mid_dB_silicon,"
                "F_min_norm_glass,F_min_norm_silicon,"
                "F_max_glass,F_max_silicon,"
                "T_coupling_741kHz\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        int n = (int)(L / 1.0e-3 + 0.5);
        if (n < 1) n = 1;
        double freq = n * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;
        double node_sp = L / (double)n;

        /* Q factor: Q ~ pi*n / (1 - r^2 + 2*alpha*L) */
        double Q_g = M_PI * n / (1.0 - r_glass*r_glass + 2.0*alpha*L);
        double Q_s = M_PI * n / (1.0 - r_silicon*r_silicon + 2.0*alpha*L);

        /* Contrast at midpoint */
        double S_g = sw_ratio(L/2.0, L, r_glass, alpha);
        double S_s = sw_ratio(L/2.0, L, r_silicon, alpha);
        double C_g = 10.0 * log10(dmax(contrast_from_swr(S_g), 1e-15));
        double C_s = 10.0 * log10(dmax(contrast_from_swr(S_s), 1e-15));

        /* Find min and max nodal force */
        double F_min_g = 1e30, F_max_g = 0;
        double F_min_s = 1e30, F_max_s = 0;
        int in;

        for (in = 0; in < n; in++) {
            double x = (in + 0.5) * node_sp;
            double Fg = fabs(radiation_force(x, L, freq, r_glass, alpha));
            double Fs = fabs(radiation_force(x, L, freq, r_silicon, alpha));
            if (Fg < F_min_g) F_min_g = Fg;
            if (Fg > F_max_g) F_max_g = Fg;
            if (Fs < F_min_s) F_min_s = Fs;
            if (Fs > F_max_s) F_max_s = Fs;
        }

        /* Coupling at 741 kHz, 10 um couplant */
        double T_coupl = coupling_T_power(741.5e3, 10e-6);

        fprintf(fp, "%.4f,%d,%.4f,%.2f,%.2f,%.4f,%.4f,"
                    "%.6f,%.6f,%.8e,%.8e,%.6f\n",
                L*1e3, n, freq/1e3,
                Q_g, Q_s, C_g, C_s,
                (F_max_g > 0) ? F_min_g/F_max_g : 0,
                (F_max_s > 0) ? F_min_s/F_max_s : 0,
                F_max_g, F_max_s,
                T_coupl);
    }

    fclose(fp);
    printf("  Written: %s (%d lengths)\n", fname, NL);
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    double r_glass, r_silicon;
    double f_target, lambda_f, lambda_2;
    int mkdir_result;

    /* Create output directory */
    mkdir_result = MKDIR(OUTPUT_DIR);
    if (mkdir_result != 0) {
        if (errno != EEXIST) {
            fprintf(stderr, "Error: cannot create %s\n", OUTPUT_DIR);
            return 1;
        }
    }

    /* Compute derived impedances */
    Z_f     = rho_f * c_f;
    Z_glass = rho_glass * c_glass;
    Z_si    = rho_si * c_si;
    Z_1     = rho_1 * c_1;
    Z_2     = rho_2 * c_2;
    Z_3     = rho_3 * c_3;

    /* BEGIN ADDITION: Compute Gor'kov contrast factor
     *
     * f1 (monopole) = 1 - kappa_p/kappa_f
     *               = 1 - (rho_f * c_f^2) / (rho_p * c_p^2)
     *
     * f2 (dipole)   = 2*(rho_p - rho_f) / (2*rho_p + rho_f)
     *
     * Phi = f1/3 + f2/2
     *
     * For typical cells in water (rho_p~1050, c_p~1550):
     *   f1 ~ 1 - (1000*1483^2)/(1050*1550^2) ~ 0.13
     *   f2 ~ 2*(1050-1000)/(2*1050+1000)     ~ 0.032
     *   Phi ~ 0.13/3 + 0.032/2               ~ 0.060
     *
     * Note: literature often cites Phi ~ 0.17 for polystyrene
     * beads (rho=1050, c=2350). For soft cells, Phi is lower
     * because c_p is closer to c_f. This matters! */
    {
        double kappa_f = 1.0 / (rho_f * c_f * c_f);   /* fluid compressibility */
        double kappa_p = 1.0 / (rho_p * c_p * c_p);   /* particle compressibility */
        f1_monopole = 1.0 - kappa_p / kappa_f;
        f2_dipole   = 2.0 * (rho_p - rho_f) / (2.0 * rho_p + rho_f);
        Phi_contrast = f1_monopole / 3.0 + f2_dipole / 2.0;
    }
    /* END ADDITION: Gor'kov contrast factor */

    r_glass  = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    f_target = c_f / (2.0 * 1.0e-3);  /* f for lambda/2 = 1 mm */
    lambda_f = c_f / f_target;
    lambda_2 = c_2 / f_target;

    printf("=====================================================\n");
    printf(" BENG207_2: Rigid Cavity Resonator Model\n");
    printf(" Lateral standing wave for organoid formation\n");
    printf("=====================================================\n");
    printf(" Fluid: water\n");
    printf("   Z_f           = %.2f MRayl\n", Z_f/1e6);
    printf("   c_f           = %.0f m/s\n", c_f);
    printf("-----------------------------------------------------\n");
    printf(" Target organoid spacing = 1 mm (CI electrode match)\n");
    printf("   f_target      = %.1f kHz\n", f_target/1e3);
    printf("   lambda_f      = %.3f mm\n", lambda_f*1e3);
    printf("   lambda/2      = %.3f mm\n", lambda_f*1e3/2.0);
    printf("-----------------------------------------------------\n");
    printf(" Wall reflection coefficients:\n");
    printf("   Glass:   Z = %.2f MRayl  |r| = %.4f\n",
           Z_glass/1e6, r_glass);
    printf("   Silicon: Z = %.2f MRayl  |r| = %.4f\n",
           Z_si/1e6, r_silicon);
    printf("-----------------------------------------------------\n");
    printf(" Vertical coupling (at %.0f kHz):\n", f_target/1e3);
    printf("   Z_LiNbO3     = %.2f MRayl\n", Z_1/1e6);
    printf("   Z_couplant   = %.2f MRayl\n", Z_2/1e6);
    printf("   Z_glass      = %.2f MRayl\n", Z_3/1e6);
    printf("   lambda_coupl = %.3f mm (d/lambda < 0.025)\n", lambda_2*1e3);
    printf("-----------------------------------------------------\n");
    printf(" Channel lengths: 5-20 mm -> 5-20 organoids\n");
    printf(" Human cochlea ~30 mm, CI = 22 electrodes\n");
    /* BEGIN ADDITION: Gor'kov summary + BENG207_3 hand-off */
    printf("-----------------------------------------------------\n");
    printf(" Gor'kov acoustic contrast factor:\n");
    printf("   rho_p  = %.0f kg/m^3 (cell density)\n", rho_p);
    printf("   c_p    = %.0f m/s (cell sound speed)\n", c_p);
    printf("   a      = %.0f um (organoid radius)\n", a_organoid*1e6);
    printf("   f1     = %.4f (monopole — compressibility contrast,\n",
           f1_monopole);
    printf("              cells less compressible than water: c_p > c_f)\n");
    printf("   f2     = %.4f (dipole — density contrast,\n", f2_dipole);
    printf("              cells slightly denser: rho_p > rho_f)\n");
    printf("   Phi    = %.4f (acoustic contrast factor)\n", Phi_contrast);
    printf("   Phi > 0 => organoids migrate to pressure NODES\n");
    printf("-----------------------------------------------------\n");
    printf(" Hand-off to BENG207_3 (stress estimation at p_ac = 500 kPa):\n");
    {
        double p_ref = 0.5e6;  /* 500 kPa reference pressure */
        double f_lat = 741.5e3;
        double f_vert = 20.0e6;
        double D_org = 2.0 * a_organoid;
        double lambda_vert = c_f / f_vert;

        /* Gorkov: sigma = 4 * Phi * k * a * E_ac
         *   where E_ac = p^2 / (4 * rho_f * c_f^2)               */
        double k_lat  = 2.0 * M_PI * f_lat / c_f;
        double k_vert = 2.0 * M_PI * f_vert / c_f;
        double E_ac   = p_ref * p_ref / (4.0 * rho_f * c_f * c_f);

        double sig_741k = 4.0 * Phi_contrast * k_lat  * a_organoid * E_ac;
        double sig_20M  = 4.0 * Phi_contrast * k_vert * a_organoid * E_ac;
        double g_corr   = 0.5;
        double sig_used = g_corr * sig_20M;

        printf("   sigma_gorkov(741 kHz)  = %.2f Pa  (D/lambda=%.2f, "
               "Gorkov valid)\n",
               sig_741k, D_org / (c_f / f_lat));
        printf("   sigma_gorkov(20 MHz)   = %.2f Pa  (raw, no g_corr)\n",
               sig_20M);
        printf("   sigma_0 -> BENG207_3   = %.2f Pa  (g_corr=%.1f for "
               "D/lambda=%.1f)\n",
               sig_used, g_corr, D_org / lambda_vert);
        printf("   [replaces old hand-waved alpha=0.15 -> 11.4 Pa]\n");
    }
    /* END ADDITION */
    printf("=====================================================\n\n");

    printf("Generating Part A outputs (lateral cavity)...\n\n");
    output_pressure_field();
    output_contrast_map();
    output_mode_structure();
    output_nodal_uniformity();

    printf("\nGenerating Part B output (vertical coupling)...\n\n");
    output_coupling_layer();

    printf("\nGenerating combined outputs...\n\n");
    output_frequency_sensitivity();
    output_design_space();

    printf("\nDone. 7 CSV files generated in %s/\n", OUTPUT_DIR);

#ifdef _WIN32
    printf("Press Enter to exit...");
    getchar();
#endif

    return 0;
}
