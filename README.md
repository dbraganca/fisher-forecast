# EFT Fisher forecast for future galaxy surveys

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

To run the project, you will need to install the requisites from the [python-integer-powers](https://github.com/dbraganca/python-integer-powers.git) GitHub repo.
We will use it to generate the matrices in order to quickly evaluate 1-loop integrals.
You need `Class` to generate the power spectra and the transfer functions.
You also need `numpy`, `matplotlib`, `getdist`, `scipy`, and `tabulate` to plot the posteriors.

### Installation

1. Go to a new directory and clone the repo
   ```
   git clone https://github.com/dbraganca/fisher-forecast.git
   ```
2. Create a new `paths.json` using the example by selecting the paths you want for the ctabs, $k$-dependent 1-loop bispectrum coefs, and the J-matrices for the power spectrum and bispectrum.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

First, we need to generate all the $k$-dependent loop coefficients.
To do so, run `make_EFT_kernels_with_bias.wl`.
This notebook generates the general bias kernels for the 1-loop power spectrum and the 1-loop bispectrum.
It will decompose each diagram kernel in so-called _ctabs_ containing the exponents of $q^{-2}$, $|k_1+q|^{-2}$, and $|k_2+q|^{-2}$, and the corresponding bias-dependent coefficient. 
The exponents of the _ctabs_ are saved in a folder that will be read to calculate the $J$-matrices.

To make computation more efficient, we decompose each bias-dependent coefficient into a combination of $A(k_1,k_2,k_3) \cdot B(\{b_i\})$, where $A$ is a vector only dependent in the $k$'s (not on any bias), and $B$ is only dependent on the biases (does not depend on any $k_i$). Permutations are also taken care of. 
The $k$-dependent $A$ coefficients are pre-computed and saved, in order to speed up to computation of the 1-loop bispectrum for each triangle. 
They require ~20 GB of memory to store. 

Second, go to the directory `python-integer-powers/1. source/` and run `python B1loop_bias.py`. 
This will generate and save the $J$-matrices required to compute loop integrals for the 1-loop bispectrum. The $J$-matrices also take ~20 GB for all the triangles.

Third, run `compute_Fisher_Surveys.wls` to compute Fisher matrices for each survey. 
This script loads the necessary functions from `FisherFunctions.wl` to calculate the Fisher matrices, but you still need to specify the values of the fiducial biases.
Computing the derivatives takes some time.
This script saves the Fisher matrices for each survey and for different specifications.

Finally, run `Full_Fisher_Plots.ipynb` to plot the predicted posteriors.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [ ] Feature 1
- [ ] Feature 2
- [ ] Feature 3
    - [ ] Nested Feature

See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>
