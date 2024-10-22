{ pkgs ? import <nixpkgs> {} }:
pkgs.clangStdenv.mkDerivation{
    # nativeBuildInputs is usually what you want -- tools you need to run
    name = "clang derivation";
    nativeBuildInputs = with pkgs.buildPackages; [ gcc cmake boost gnuplot];
}
