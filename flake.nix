{
  description = "Single-file linear algebra math library for C++23.";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          system = system;
        };
      in
      {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            cmake
            ninja
            clang-tools
            lldb
            pkg-config
          ];
        };

        packages.default = pkgs.stdenv.mkDerivation {
          pname = "smath";
          version = "master";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            cmake
            gtest
            copyPkgconfigItems
          ];

          pkgconfigItems = [
            (pkgs.makePkgconfigItem rec {
              name = "smath";
              version = "1";
              cflags = [ "-I${variables.includedir}" ];
              variables = rec {
                prefix = "${placeholder "out"}";
                includedir = "${prefix}/include";
              };
              description = "Single-file linear algebra math library for C++23.";
            })
          ];

          installPhase = ''
            runHook preInstall
            mkdir -p $out/include
            cp ../include/smath.hpp $out/include/
            runHook postInstall
          '';

          dontBuild = false;
          doCheck = true;

          meta = with pkgs.lib; {
            description = "Single-file linear algebra math library for C++23.";
            homepage = "https://github.com/slendidev/smath";
            license = licenses.asl20;
            platforms = platforms.all;
            maintainers = [
              {
                name = "Slendi";
                email = "slendi@socopon.com";
                github = "slendidev";
                githubId = 32436619;
              }
            ];
          };
        };
      }
    );
}
