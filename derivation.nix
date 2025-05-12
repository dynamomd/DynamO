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

  buildInputs = [
    pkgs.cmake
    pkgs.git
    pkgs.ninja
    pkgs.gcc
    pkgs.pkg-config
    pkgs.bzip2.dev
    ((pks: pks.python3.withPackages (ps: with ps; [
      numpy
    ])) pkgs)
  ];

  configurePhase = ''
    cmake .
  '';

  buildPhase = ''
    cmake --build .
  '';

  installPhase = ''
    mkdir -p $out/bin
    mv chord $out/bin
  '';
}
