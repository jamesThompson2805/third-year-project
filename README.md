# third-year-project

The final report is also included in this repository under project.pdf .

For usage of this project, either consider installation of Nix to handle management of required libraries for you or open flake.nix and read requirements from the list defined within ```packages = with pkgs; [```.

Some method for compiling cpp (gcc, clang or other) is needed alongside cmake.
The boost library needs to be accessible to the project by adjusting CMakeLists.txt, installing system wide or using Nix.
Gnuplot is needed for graph rendering and results will attempt to be viewed with the Floorp browser if PDF is selected (the pdf will be saved regardless).
Directories used in code are in style of Linux so will need to be adjusted for use with windows (replace / with \).

The directory _lib_ holds all the implemented library code for parsing, cleaning, approximation, indexing and similarity Search.
The directory _src_ code for the main, demo, evaluations and plotting.


### Abstract:
As more data is collected in the world, databases of Time Series Data become larger and increase the strain on algorithms operating on the data. Compression techniques have been developed that maintain structural properties of the data to allow large datasets to be stored and processed efficiently. Adaptive Piecewise Linear Approximation is one compression technique that approximates continuous segments of data. However computing the optimal segments is challenging and multiple algorithms across multiple disciplines have been created.

Similarity Search is a fundamental problem in Time Series Analysis, for example used to detect heart beats and find stock market trends. This codebase supports a report providing a review and empirical comparison of current Adaptive PLA algorithms. A novel approach to extend Adaptive PLA to solve Similarity Search is implemented, capable of using high compression techniques and offering unique advantages over other indexing schemes. 


