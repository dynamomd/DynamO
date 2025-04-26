{
  description = "A Nix-flake-based C/C++ development environment";

  inputs.nixpkgs.url = "https://flakehub.com/f/NixOS/nixpkgs/0.1.*.tar.gz";

  outputs = { self, nixpkgs }:
    let
      supportedSystems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forEachSupportedSystem = f: nixpkgs.lib.genAttrs supportedSystems (system: f {
        pkgs = import nixpkgs { inherit system; };
      });
      extraOutputsToInstall = [ "dev" ];
    in
    {
      devShells = forEachSupportedSystem ({ pkgs }: {
        default = pkgs.mkShell.override
          {
            # Override stdenv in order to change compiler:
            #stdenv = pkgs.clangStdenv;
          }
          {
            packages = with pkgs; [
              clang
              clang-tools
              cmake
              ninja
              doxygen
              vcpkg
              vcpkg-tool
              boost.dev
              eigen
              bzip2.dev
              python3
              gtkmm3.dev
              ffmpeg.dev
              freeglut
              glew
              cairomm
              pkg-config
            ] ++ (if system == "aarch64-darwin" then [ ] else [ gdb ]);
          };
      });
    };
}

