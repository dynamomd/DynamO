{
  description = "PyDynamO / DynamO build environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        pydynamo = pkgs.callPackage ./derivation.nix {};
      in
      {
        packages.default = pydynamo;
        packages.pydynamo = pydynamo;

        devShells.default = pkgs.mkShell {
          inputsFrom = [ pydynamo ];
        };
      }
    );
}
