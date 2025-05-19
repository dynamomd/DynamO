# see https://unix.stackexchange.com/questions/717168/how-to-package-my-software-in-nix-or-write-my-own-package-derivation-for-nixpkgs
# To install `nix-env -u -f default.nix`
# To develop `nix-shell` (will build the shell with dependencies)
# To test build `nix-build`
{ pkgs, python3 }:
python3.pkgs.buildPythonPackage rec {
  name = "pydynamo";
  src = ./.;
  pyproject = true;

  dontUseCmakeConfigure = true;
  
  propagatedBuildInputs = with python3.pkgs; [
      scikit-build-core
      numpy
      alive-progress
      uncertainties
      pandas
      scipy
      freud
      networkx
      pytest
  ];

  nativeBuildInputs = with pkgs; [
    cmake
    ninja
    git
    gcc
    pkg-config
    clang-tools
    wrapGAppsHook3
  ]  ++ propagatedBuildInputs;
  
  buildInputs = with pkgs; [
    # Basic build dependencies
    bzip2.dev
    boost.dev
    eigen
    
    # Visualiser
    libGL.dev
    gtkmm3.dev
    ffmpeg.dev
    freeglut
    glew
    cairomm.dev
    libpng
  ];
}
