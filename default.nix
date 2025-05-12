{
  pkgs ? import <nixpkgs> {}
}:
pkgs.stdenv.mkDerivation rec {
  pname = "dynamomd";
  version = "1.7.0";
  #src = pkgs.fetchgit {
  #  url = "https://github.com/dynamomd/dynamo";
  #  rev =  "f1fd582cd9ed09b52b7c99cdbf00cd1cebb3958c"; # "refs/heads/nix";
  #  sha256 = "sha256-S3D15QD4NTVSM6PR7xqRQj7yvpq2MQj4WmHzspDKzTI=";
  #};

  src = ./.;
  
  buildInputs = with pkgs; [
    # Basic build dependencies
    cmake
    git
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

  #configurePhase = ''
  #  cmake .
  #'';
  #
  #buildPhase = ''
  #  cmake --build . -j32
  #'';
  #
  #installPhase = ''
  #  cmake --install .
  #'';
}
