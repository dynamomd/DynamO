{ pkgs ? import <nixpkgs> {
  config.allowUnfree = true;
} }:
let
  envname = "platformio-fhs";
in
(pkgs.buildFHSUserEnv {
  name = envname;
  targetPkgs = pkgs: (with pkgs; [
    bashInteractive # Needed to fix shell in vscode
    git
    ((pks: pks.python3.withPackages (ps: with ps; [
      numpy
    ])) pkgs)

    # Basic build dependencies
    cmake
    
    # VSCode with all the extensions I use
    (vscode-with-extensions.override {
      vscodeExtensions = with pkgs.vscode-extensions; [
      tuttieee.emacs-mcx
	    ms-vscode.cpptools
	    ms-vscode.cpptools-extension-pack
	    xaver.clang-format
	    ms-vscode.cmake-tools
	    #brobeson.ctest-lab
	    batisteo.vscode-django
	    grapecity.gc-excelviewer
	    github.vscode-github-actions
	    github.copilot
	    ms-vscode.live-server
	    #elagil.pre-commit-helper
	    ms-python.python
	    ms-python.debugpy
      ms-vscode-remote.remote-ssh
      dbaeumer.vscode-eslint
      tomoki1207.pdf
      ];
    })
  ]);
}).env
