# CV Optimization via Kinetic Rate Minimization

This project aims to optimize collective variables (CVs) by minimizing the kinetic rate. The provided code, written in Fortran, utilizes a Monte Carlo Metropolis optimization algorithm to generate an optimized set of weights for a linear combination of the CVs. The optimized weights are obtained by iteratively perturbing an initial random combination and accepting perturbations that minimize the kinetic rate.

## Theory

The kinetic rate (k_{AB}) is calculated using the following formula:

k_{AB}^{-1} = \int_{q_0}^{b} dx \frac{e^{\beta F(x)}}{D(x)} \int_{a}^{x} dy e^{-\beta F(y)}

Where:
- F(x) represents the free energy.
- D(x) represents the diffusion.
- The integral limits are defined as follows:
  - q_0: Lower limit for the integration.
  - b: Upper limit for the integration.
  - a: Lower limit for the inner integration.

The values for F(x) and D(x) are obtained by training an overdamped Langevin model in the high friction regime via likelihood maximization. For more details about the training process and methodology, refer to the publication: "Free Energy Landscapes, Diffusion Coefficients, and Kinetic Rates from Transition Paths" published in J. Chem. Theory Comput. 2022, 18, 8, 4639–4648.

## Getting Started

To use this code, follow the steps below:

### Prerequisites

Make sure you have the following dependencies installed on your system:

- Fortran compiler (e.g., GNU Fortran)
- optLE executable (see below)

### OptLE Code

This code relies on the optLE code, which needs to be downloaded from the [optLE GitHub repository](https://github.com/physix-repo/optLE). Please refer to their README file for instructions on how to obtain and compile the optLE executable.


### Compilation

Compile the `optCol.f90` code using a Fortran compiler. For example, using the GNU Fortran compiler:

gfortran optCol.f90 -o optCol.x


### Usage

1. Prepare your input file:
   - Create a file named `colvar_x1` containing the input set of CVs, with each CV represented as a column and time as the first column

2. Run the compiled executable:
./optCol.x 


3. Review the output:
- The optimized weights for the linear combination of CVs will be stored in the first columns of the output file `tau`.
- The last column of the `tau` file contains the mean first passage time, which is the inverse of the rate.

## Contributing

Contributions to this project are welcome! If you have any ideas, suggestions, or improvements, please feel free to open an issue or submit a pull request.

## Acknowledgments

- More details can be found in our arXiv preprint: arXiv preprint arXiv:2302.12497
- This work is done by Line Mouaffac, Karen Palacio-Rodriguez and Fabio Pietrucci

## Contact

If you have any questions or need further assistance, please feel free to contact Line Mouaffac at line.mouaffac@sorbonne-universite.fr .




