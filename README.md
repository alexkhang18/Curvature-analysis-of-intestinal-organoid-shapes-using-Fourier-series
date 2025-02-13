# Curvature-analysis-of-intestinal-organoid-shapes-using-Fourier-series
This code takes in an intestinal organoid segmentation and constructs a Fourier series representation of the boundary. Curvature along the boundary representation is then computed and interpolated within the organoid body to attribute curvature values to regions within the organoid. There are three files included that are required for the code to run: 'organoid_brightfield.tif', 'organoid_segmentation.png', and 'lysozyme_segmentation.png'. Ensure all the necessary files are in the same directory as the main code file titled 'curvature_2D_fourier_paneth_interpolation.m' and the functions 'fourier_series_fit.m' and 'fourier_series.m'. Allow a few minutes for the code to run and generate the following images:

![Picture1](https://github.com/user-attachments/assets/9fb3472c-21c9-4448-9e13-316094b11975)

The code should also generate the plots below which show the arc-length of the boundary of the organoid that falls into a specific curvature range and the probability density estimate of the curvature values associated with lysozyme positive pixels:

![Picture2](https://github.com/user-attachments/assets/bb425555-76d7-4472-b7ee-0e6b785327ad)

The code relies on the following segmentations of the organoid body where blue marks the boundary: 

![organoid_brightfield_segmentation_check](https://github.com/user-attachments/assets/7b3871d5-2fe0-440e-8a68-57eac4507e51)

The code also relies on the following segmentation of lysozyme positive pixels (note that the positive pixels of the necrotic core are ignored): 

![lysozyme_segmentation_check](https://github.com/user-attachments/assets/549e24a8-9b7a-4140-b37b-f1649c45dd92)
