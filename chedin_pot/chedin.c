#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

/* 
   Helper to apply minimum image convention for an orthorhombic cell.
   Brings the distance component `dx` back into the [-L/2, L/2] range.
*/
inline double wrap_pbc(double dx, double L) {
    if (L <= 0.0) return dx; // No PBC along this axis if length is 0
    return dx - L * round(dx / L);
}

/*
   Core computational loop with Periodic Boundary Conditions.
*/
void compute_chedin_pbc(int n_atoms, double* pos, double* cell, double* energy, double* f) {

    double pi = 3.14159265358979323846;
    double pot = 0.0;

    // Extract box lengths for X, Y, and Z axes (assuming diagonal/orthorhombic cell)
    double Lx = cell[0 * 3 + 0]; // cell[0][0]
    double Ly = cell[1 * 3 + 1]; // cell[1][1]
    double Lz = cell[2 * 3 + 2]; // cell[2][2]

    // Initialize forces array to zero
    for (int i = 0; i < n_atoms * 3; i++) f[i] = 0.0;

    // Cartesian distances
    double dx21 = pos[3] - pos[0];
    double dy21 = pos[4] - pos[1];
    double dz21 = pos[5] - pos[2];

    double dx23 = pos[6] - pos[0];
    double dy23 = pos[7] - pos[1];
    double dz23 = pos[8] - pos[2];

    // Minimum Image Convention (PBC wrap)
    dx21 = wrap_pbc(dx21, Lx);
    dy21 = wrap_pbc(dy21, Ly);
    dz21 = wrap_pbc(dz21, Lz);

    dx23 = wrap_pbc(dx23, Lx);
    dy23 = wrap_pbc(dy23, Ly);
    dz23 = wrap_pbc(dz23, Lz);

    // Distances
    double r21 = sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21);
    double r23 = sqrt(dx23*dx23 + dy23*dy23 + dz23*dz23);

    // Cross product and bond angle
    double ncross_x = dy21 * dz23 - dz21 * dy23;
    double ncross_y = dz21 * dx23 - dx21 * dz23;
    double ncross_z = dx21 * dy23 - dy21 * dx23;

    double nabs = sqrt(ncross_x * ncross_x + ncross_y * ncross_y + ncross_z * ncross_z);
    double alpha = asin(nabs/(r21*r23));

    // Equilibrium C-O distance
    double re = 1.161;
    // Quasi-normal coordinates parameters
    double a1 = 20.294; 
    double a2 = 7.665; 
    double a3 = 14.248;

    // Dimensionless coordinates
    double xi1 = 1 - re/r21;
    double xi2 = alpha;
    double xi3 = 1 - re/r23;

    // Quasi-normal symmetry coordinates (Chedin)
    double rq1 = a1*(xi1+xi3);
    double rq2 = a2*xi2;
    double rq3 = a3*(xi1-xi3);

    // Quadratic terms in potential, 3 parameters
    double par2[3] = {8.49441, 4.01316, 14.7322};
    pot += par2[0]*rq1*rq1 + par2[1]*rq2*rq2 + par2[2]*rq3*rq3;

    // Cubic & 4th order terms in potential, 3 + 6 = 9 parameters
    double par34[9] = {-0.169658, -0.349882, -0.915698, -0.00309579, 0.00762742, -0.0349418, 0.00988316, -0.0332545, -0.0148403};
    pot += par34[0]*rq1*rq1*rq1 + par34[1]*rq1*rq2*rq2 + par34[2]*rq1*rq3*rq3
         + par34[3]*rq1*rq1*rq1*rq1 + par34[4]*rq1*rq1*rq2*rq2 
         + par34[5]*rq1*rq1*rq3*rq3 + par34[6]*rq2*rq2*rq2*rq2 
         + par34[7]*rq2*rq2*rq3*rq3 + par34[8]*rq3*rq3*rq3*rq3;

    // 5th order terms in potential, 6 parameters
    double par5[6] = {0.000567461, 0.000109297, 0.00727883, -0.00135831, -0.00149995, 0.00120956};
    pot += par5[0]*rq1*rq1*rq1*rq1*rq1 + par5[1]*rq1*rq1*rq1*rq2*rq2 
         + par5[2]*rq1*rq1*rq1*rq3*rq3 + par5[3]*rq1*rq2*rq2*rq2*rq2 
         + par5[4]*rq1*rq2*rq2*rq3*rq3 + par5[5]*rq1*rq3*rq3*rq3*rq3;

    // 6th order terms in potential, 10 parameters
    double par6[10] = {2.53031e-5, -1.00491e-3, -1.23064e-3, -1.06667e-4, -5.10572e-4, -1.67567e-3, 8.90746e-5, -1.3157e-4, -1.08281e-5, -1.38092e-4};
    pot += par6[0]*rq1*rq1*rq1*rq1*rq1*rq1 + par6[1]*rq1*rq1*rq1*rq1*rq2*rq2 
         + par6[2]*rq1*rq1*rq1*rq1*rq3*rq3 + par6[3]*rq1*rq1*rq2*rq2*rq2*rq2 
         + par6[4]*rq1*rq1*rq2*rq2*rq3*rq3 + par6[5]*rq1*rq1*rq3*rq3*rq3*rq3
         + par6[6]*rq2*rq2*rq2*rq2*rq2*rq2 + par6[7]*rq2*rq2*rq2*rq2*rq3*rq3
         + par6[8]*rq2*rq2*rq3*rq3*rq3*rq3 + par6[9]*rq3*rq3*rq3*rq3*rq3*rq3;


    // Derivatives over quasi-normal coordinates
    double dudrq1 = 2*par2[0]*rq1
             + 3*par34[0]*rq1*rq1 + 4*par34[3]*rq1*rq1*rq1 
               + (par34[1] + 2*par34[4]*rq1)*rq2*rq2 + (par34[2] + 2*par34[5]*rq1)*rq3*rq3
             + 5*par5[0]*rq1*rq1*rq1*rq1 + 3*par5[1]*rq1*rq1*rq2*rq2 + par5[3]*rq2*rq2*rq2*rq2 
  	           + (3*par5[2]*rq1*rq1 + par5[4]*rq2*rq2)*rq3*rq3 + par5[5]*rq3*rq3*rq3*rq3
             + 2*rq1*(3*par6[0]*rq1*rq1*rq1*rq1 + 2*par6[1]*rq1*rq1*rq2*rq2 
  	           + par6[3]*rq2*rq2*rq2*rq2 + (2*par6[2]*rq1*rq1 + par6[4]*rq2*rq2)*rq3*rq3
  	           + par6[5]*rq3*rq3*rq3*rq3);
    double dudrq2 = 2*par2[1]*rq2
             + 2*rq2*(rq1*(par34[1] + par34[4]*rq1) + 2*par34[6]*rq2*rq2 + par34[7]*rq3*rq3)
             + 2*rq1*rq2*(par5[1]*rq1*rq1 + 2*par5[3]*rq2*rq2 + par5[4]*rq3*rq3)
             + 2*rq2*(par6[1]*rq1*rq1*rq1*rq1 + 2*par6[3]*rq1*rq1*rq2*rq2 
  	           + 3*par6[6]*rq2*rq2*rq2*rq2 + (par6[4]*rq1*rq1 + 2*par6[7]*rq2*rq2)*rq3*rq3
  	           + par6[8]*rq3*rq3*rq3*rq3);
    double dudrq3 = 2*par2[2]*rq3
             + 2*rq3*(rq1*(par34[2] + par34[5]*rq1) + par34[7]*rq2*rq2 + 2*par34[8]*rq3*rq3)
             + 2*rq1*rq3*(par5[2]*rq1*rq1 + par5[4]*rq2*rq2 + 2*par5[5]*rq3*rq3)
             + 2*rq3*(par6[2]*rq1*rq1*rq1*rq1 + par6[4]*rq1*rq1*rq2*rq2 
  	           + par6[6]*rq2*rq2*rq2*rq2 + 2*(par6[5]*rq1*rq1 + par6[8]*rq2*rq2)*rq3*rq3
  	           + 3*par6[9]*rq3*rq3*rq3*rq3);

    // Forces

    f[3] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*dx21;
    f[4] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*dy21;
    f[5] = -(dudrq1*a1+dudrq3*a3)*re/(r21*r21*r21)*dz21;
    
    f[6] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*dx23;
    f[7] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*dy23;
    f[8] = -(dudrq1*a1-dudrq3*a3)*re/(r23*r23*r23)*dz23;

    if (nabs > 0) {
        f[3] -= dudrq2*a2/(nabs/r23)*(dx23/r23-dx21*cos(pi-alpha)/r21);
        f[4] -= dudrq2*a2/(nabs/r23)*(dy23/r23-dy21*cos(pi-alpha)/r21);
        f[5] -= dudrq2*a2/(nabs/r23)*(dz23/r23-dz21*cos(pi-alpha)/r21);
        
        f[6] -= dudrq2*a2/(nabs/r21)*(dx21/r21-dx23*cos(pi-alpha)/r23);
        f[7] -= dudrq2*a2/(nabs/r21)*(dy21/r21-dy23*cos(pi-alpha)/r23);
        f[8] -= dudrq2*a2/(nabs/r21)*(dz21/r21-dz23*cos(pi-alpha)/r23);
    }

    f[0] = - f[3] - f[6];
    f[1] = - f[4] - f[7];
    f[2] = - f[5] - f[8];

    *energy = pot;
}


/* Python Interface Wrapper */
static PyObject* method_compute(PyObject* self, PyObject* args) {
    PyArrayObject *pos_array = NULL;
    PyArrayObject *cell_array = NULL;
    PyArrayObject *forces_array = NULL;

    // Expecting 3 array arguments from Python: positions, cell matrix, forces
    if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &pos_array, 
                                          &PyArray_Type, &cell_array, 
                                          &PyArray_Type, &forces_array)) {
        return NULL;
    }

    int n_atoms = (int)PyArray_DIM(pos_array, 0);

    // Extract raw pointers
    double *pos = (double*)PyArray_DATA(pos_array);
    double *cell = (double*)PyArray_DATA(cell_array);
    double *forces = (double*)PyArray_DATA(forces_array);
    double energy = 0.0;

    // Compute potential
    compute_chedin_pbc(n_atoms, pos, cell, &energy, forces);

    return PyFloat_FromDouble(energy);
}

// Module registration boilerplates
static PyMethodDef MyMethods[] = {
    {"compute", method_compute, METH_VARARGS, "Calculates energy and forces with PBC."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mymodule = {
    PyModuleDef_HEAD_INIT, "ChedinPot_calculator", NULL, -1, MyMethods
};

PyMODINIT_FUNC PyInit_ChedinPot_calculator(void) {
    import_array(); 
    return PyModule_Create(&mymodule);
}
