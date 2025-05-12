{
  pkgs ? import <nixpkgs> {}
}:
pkgs.stdenv.mkDerivation rec {
  pname = "dynamomd";
  version = "1.6.0";
  src = pkgs.fetchgit {
    url = "https://github.com/dynamomd/dynamo";
    rev = "refs/heads/nix";
    sha256 = "sha256-p8+OoNBW7VABgIWXme6iLWEkiPe7v9yZqPveN4A+hKY=";
  };

  buildInputs = with pkgs; [
    # Basic build dependencies
    cmake
    git
    ninja
    gcc
    pkg-config
    bzip2.dev
    boost.dev
    clang-tools
    eigen
    ((pks: pks.python3.withPackages (ps: with ps; [
      numpy
      #alive-progress
      uncertainties
      pandas
      scipy
      freud
      networkx
    ])) pkgs)

    # Visualiser
    gtkmm3.dev
    ffmpeg.dev
    freeglut
    glew
    cairomm
    libpng
  ];

  configurePhase = ''
    cmake .
  '';

  buildPhase = ''
    cmake --build .
  '';

  installPhase = ''
    cmake --install .
  '';
}
