#!/usr/bin/env bash
# Copyright (2026) Yuanmai Zheng
# SPDX-License-Identifier: GPL-2.0-only
set -euo pipefail
export LC_ALL=C
package_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
[[ $# -eq 1 ]] || { echo 'Usage: install.sh TARGET_CHECKOUT' >&2; exit 2; }
target_root=$(cd -- "$1" && pwd -P)
digest() { sha256sum -- "$1" | cut -d ' ' -f1; }
fail() { printf '%s\n' "$1" >&2; exit 2; }
[[ "$target_root" != "$package_root" ]] || fail 'SPINGNN_INSTALL_INVALID_TARGET'
installed=1
while read -r expected relative; do
  if [[ ! -f "$target_root/$relative" || -L "$target_root/$relative" ]] || [[ $(digest "$target_root/$relative") != "$expected" ]]; then installed=0; break; fi
done <<'FINAL_SHA256'
732bebce98645ad170a574834d307b3dc234ac523d94e9d57a297a99cc868e15 lammps/NOTICE.md
444bbcdf5dcd12fa5becf8e53811f2e2a068a70ef064be6c001bbff421799cb0 lammps/cmake/CMakeLists.txt
207221363ed61e601a96636fba68a4b7c3c99b981753def4d9e857c3ebb49145 lammps/cmake/Modules/Packages/SPINGNN.cmake
3f9d66844a03f0b97959e9be770e535971076e77999c0a4dec45085b686e6d1a lammps/doc/src/SPINGNN.rst
55bafd18d43ad742eaf80c3ccab2b7485f6ca50a7fcac2f487cb943246d425ef lammps/doc/src/fix_spin_lpsib.rst
d26b7242dfa468a361fb00606cbcfa3719dc5b67e17c8104cf4247c0f7ae4fd4 lammps/doc/src/fix_spin_sib.rst
f2603643e2f377893f979042a8d60b3329769c88388c83aa4f4579904d6a7719 lammps/doc/src/fix_spin_sib_npt.rst
8894d0db2c7be83e6558477cd092115c868148f76f119a1145e76502f189a433 lammps/doc/src/pair_spin_allegro_spin.rst
e2343cbdc0d91d81396bb8e1c6f2bca234f956ed11fce165adbf3866c60a01af lammps/doc/src/pair_spin_nequip_spin.rst
b8b3cceeacf32a3e0b26e5db223e1f310901244d2334a160146889e749c07187 lammps/examples/SPINGNN/README-C5.md
dba052daeff0d6b2ce61ff48a31a1c1ab42fb091221036642354dfb2fad254eb lammps/examples/SPINGNN/in.crsbr-lpsib
a91afa3276baa7956c675d756049896c659489ac4fd85b4a88aa49306fac6a8b lammps/examples/SPINGNN/in.fe-lpsib-centered
4cce5b13fb4b41e6426838bf36bb65b56762ca7c067720590ea0308d1034ac5d lammps/examples/SPINGNN/in.fe-lpsib-frozen
3b67241b833ebec72432ed7c10dff7843b1858068fae8008cea352f1d20218ed lammps/examples/SPINGNN/in.fe-lpsib-moving-nve
0beedcd7ea4207d0d60344f003455c0c31e7c863c5e06a1ff557c0f0e084f016 lammps/examples/SPINGNN/in.fe-lpsib-moving-nvt
3ac8df2322707c282d5db737881f58243397e0d6cfe31053e1979e2f6a52ec5a lammps/examples/SPINGNN/in.fe-lpsib-npt
ef62417cc80de9e17131b657fc7ad01b37304574ef8b98796550ba384bf01403 lammps/examples/SPINGNN/in.fe-sib-npt
6716ea449f24d95cc6c0b42742c48ebb02d9fd4a7cceaf42feba7d47da27ba78 lammps/examples/SPINGNN/in.fe-sib-nve
50d99a394d1e66cd268acef048a5bb8a6006659e4f21b5d65a52fd68caad74f8 lammps/examples/SPINGNN/in.fe-sib-nvt
ea732af643dd333b0456ad6763ce5f239ee40b9e73af6731fbda5f7a2d504676 lammps/src/Depend.sh
e6060ec0bf8d6adb4ba5f94aaa7b8f00b07169f8ed70eed8a580143d606391f8 lammps/src/Makefile
e72616b27d98bc01733d71ce7cb6acc0a40818acb06fc961aaca4eb40cfd7e33 lammps/src/SPIN/fix_spin_sib.cpp
bab7601c3c48b5ded1036d4df4b9f5d828e10d988ba7b1a5f8bafeb8ae269165 lammps/src/SPIN/pair_spin.cpp
aae06f610f92bece7742797fd065121741a71830869fd1f67a274c1fc6efe324 lammps/src/SPINGNN/CHANGELOG.md
bf67e41443cb5bf024d4113051f58f26f9c42ab706a2138ce3710d7fea367f3a lammps/src/SPINGNN/Install.sh
67754898a77a45ffec5e3934e6b5dcba8049b2aaf2cb215c4609af793ca10aad lammps/src/SPINGNN/KNOWN-LIMITATIONS.md
12df52ee7e77c83f6456022cc2b27e767126500d081ab85ac28bbe5843cab554 lammps/src/SPINGNN/README
0c125adbf8e9bd32062405a81074d77496ed58cfaaabafcaef94da0f8c161319 lammps/src/SPINGNN/doc-anchors.json
abaf6c22e9d5e2e900aebc4555946e67618e12cb59bf0670ca11954756d3a1e7 lammps/src/SPINGNN/fix_spin_lpsib.cpp
4d78270e1e491e5485f86646fa802bbb2849e06a7f8c932529427842a7a7f38f lammps/src/SPINGNN/fix_spin_lpsib.h
3f4aa405d9e19f94e22479ee545c4b22224db203c68ab35b83c2362d8960ca12 lammps/src/SPINGNN/pair_allegro_spin.cpp
509942afa6643a769287e14dfd37fe95c3dc26c35d56c09a6ed67e82f309a77e lammps/src/SPINGNN/pair_allegro_spin.h
29e41c2307c5e2cb7b36e715f3e0578d1c0c5fe63e5879c878691bf79e3c1dcf lammps/src/SPINGNN/pair_nequip_spin.cpp
330a69e10c224f36c1e43c6f83388ddc562e6380efe9bfbf9b19ac39e5d52ae9 lammps/src/SPINGNN/pair_nequip_spin.h
173e52ec260cc5de6d0d291b1aebea6a2733ed868d3f92a5f314d25015c9b3a8 lammps/src/fix_nh.cpp
037bfb80025a8686bf4766d16ca5f2afbccdc72a84d00ffa42c3a8bcfff93add lammps/src/observe85.h
f6790805b27da70d642e6c26fc33b6cd05b9856066c00dea42e46209d17c05b1 lammps/src/pair_hybrid.h
7caf34ceeb1ae32f2ce908259d2e35b36ec81fa4b7f4a72f457c2468973702b8 lammps/src/set.cpp
65e55ed653c24f86ed7ad3078c00f667672f246436f721b0b2f9209ba409f2d0 lammps/src/verlet.cpp
FINAL_SHA256
[[ ! -e "$target_root/lammps/src/SPIN/pair_allegro_spin.cpp" && ! -L "$target_root/lammps/src/SPIN/pair_allegro_spin.cpp" ]] || installed=0
[[ ! -e "$target_root/lammps/src/SPIN/pair_allegro_spin.h" && ! -L "$target_root/lammps/src/SPIN/pair_allegro_spin.h" ]] || installed=0
[[ ! -e "$target_root/lammps/src/SPIN/pair_nequip_spin.cpp" && ! -L "$target_root/lammps/src/SPIN/pair_nequip_spin.cpp" ]] || installed=0
[[ ! -e "$target_root/lammps/src/SPIN/pair_nequip_spin.h" && ! -L "$target_root/lammps/src/SPIN/pair_nequip_spin.h" ]] || installed=0
[[ $installed -eq 0 ]] || fail 'SPINGNN_INSTALL_ALREADY_INSTALLED'
while read -r expected relative; do
  [[ -f "$target_root/$relative" && ! -L "$target_root/$relative" ]] || fail "SPINGNN_INSTALL_BASE_MISMATCH:$relative"
  [[ $(digest "$target_root/$relative") == "$expected" ]] || fail "SPINGNN_INSTALL_BASE_MISMATCH:$relative"
done <<'BASE_SHA256'
9f8bd17db5982d05d18fca77a61a529a2ac5d68447d5349ff7904292cb18b744 lammps/NOTICE.md
00b4bfd35d879d35e2d5ca3468bb9f5e4fe4c1ee6a4ad26f382c481b57764613 lammps/cmake/CMakeLists.txt
787190d2251c538a588d0c2db048f6b214db1a3fc224f623797fa3cfeb7954ae lammps/src/Depend.sh
b123903c95a6c2caded11a97ab0561f84f3ad79b1643aa5ac235c0e2fb5b3166 lammps/src/Makefile
59effdf8914f23a502048383ca36d60e0074df2e02475a7c1cab9f3034b11101 lammps/src/SPIN/fix_spin_sib.cpp
3f4aa405d9e19f94e22479ee545c4b22224db203c68ab35b83c2362d8960ca12 lammps/src/SPIN/pair_allegro_spin.cpp
509942afa6643a769287e14dfd37fe95c3dc26c35d56c09a6ed67e82f309a77e lammps/src/SPIN/pair_allegro_spin.h
29e41c2307c5e2cb7b36e715f3e0578d1c0c5fe63e5879c878691bf79e3c1dcf lammps/src/SPIN/pair_nequip_spin.cpp
330a69e10c224f36c1e43c6f83388ddc562e6380efe9bfbf9b19ac39e5d52ae9 lammps/src/SPIN/pair_nequip_spin.h
805780edf90b8496eb8cc9842274199b67628b0a135b481280a88b6a03e1e311 lammps/src/SPIN/pair_spin.cpp
6e720eef78454c3fbd48c2d7984fe295afcbb780cc30cba90634def78669b318 lammps/src/fix_nh.cpp
5a4933356f5a875e5ed9737596ebb86daf207e51600435cd2d4986e8dce021ed lammps/src/pair_hybrid.h
d478753cbf27987086548560a8a1b0d6eef925bdd0fad3d80a0dab967b4d9e46 lammps/src/set.cpp
93534cd779afc247e2cbd142a3c38c215299fa52c428da5e84d80c0b176ea5d2 lammps/src/verlet.cpp
BASE_SHA256
while read -r relative; do
  [[ ! -e "$target_root/$relative" && ! -L "$target_root/$relative" ]] || fail "SPINGNN_INSTALL_TARGET_EXISTS:$relative"
done <<'DESTINATIONS'
lammps/cmake/Modules/Packages/SPINGNN.cmake
lammps/doc/src/SPINGNN.rst
lammps/doc/src/fix_spin_lpsib.rst
lammps/doc/src/fix_spin_sib.rst
lammps/doc/src/fix_spin_sib_npt.rst
lammps/doc/src/pair_spin_allegro_spin.rst
lammps/doc/src/pair_spin_nequip_spin.rst
lammps/examples/SPINGNN/README-C5.md
lammps/examples/SPINGNN/in.crsbr-lpsib
lammps/examples/SPINGNN/in.fe-lpsib-centered
lammps/examples/SPINGNN/in.fe-lpsib-frozen
lammps/examples/SPINGNN/in.fe-lpsib-moving-nve
lammps/examples/SPINGNN/in.fe-lpsib-moving-nvt
lammps/examples/SPINGNN/in.fe-lpsib-npt
lammps/examples/SPINGNN/in.fe-sib-npt
lammps/examples/SPINGNN/in.fe-sib-nve
lammps/examples/SPINGNN/in.fe-sib-nvt
lammps/src/SPINGNN/CHANGELOG.md
lammps/src/SPINGNN/Install.sh
lammps/src/SPINGNN/KNOWN-LIMITATIONS.md
lammps/src/SPINGNN/README
lammps/src/SPINGNN/doc-anchors.json
lammps/src/SPINGNN/fix_spin_lpsib.cpp
lammps/src/SPINGNN/fix_spin_lpsib.h
lammps/src/SPINGNN/pair_allegro_spin.cpp
lammps/src/SPINGNN/pair_allegro_spin.h
lammps/src/SPINGNN/pair_nequip_spin.cpp
lammps/src/SPINGNN/pair_nequip_spin.h
lammps/src/observe85.h
DESTINATIONS
while read -r relative; do
  parent="$target_root/$relative"
  while [[ "$parent" != "$target_root" ]]; do
    [[ ! -L "$parent" ]] || fail "SPINGNN_INSTALL_SYMLINK_PARENT:$relative"
    [[ ! -e "$parent" || -d "$parent" ]] || fail "SPINGNN_INSTALL_INVALID_PARENT:$relative"
    parent=$(dirname -- "$parent")
  done
done <<'PARENTS'
lammps
lammps/cmake
lammps/cmake/Modules
lammps/cmake/Modules/Packages
lammps/doc
lammps/doc/src
lammps/examples
lammps/examples/SPINGNN
lammps/src
lammps/src/SPIN
lammps/src/SPINGNN
PARENTS
while read -r expected relative; do
  [[ -f "$package_root/$relative" && ! -L "$package_root/$relative" ]] || fail "SPINGNN_INSTALL_PAYLOAD_MISMATCH:$relative"
  [[ $(digest "$package_root/$relative") == "$expected" ]] || fail "SPINGNN_INSTALL_PAYLOAD_MISMATCH:$relative"
done <<'PAYLOAD_SHA256'
207221363ed61e601a96636fba68a4b7c3c99b981753def4d9e857c3ebb49145 cmake/SPINGNN.cmake
3f9d66844a03f0b97959e9be770e535971076e77999c0a4dec45085b686e6d1a doc/SPINGNN.rst
55bafd18d43ad742eaf80c3ccab2b7485f6ca50a7fcac2f487cb943246d425ef doc/fix_spin_lpsib.rst
d26b7242dfa468a361fb00606cbcfa3719dc5b67e17c8104cf4247c0f7ae4fd4 doc/fix_spin_sib.rst
f2603643e2f377893f979042a8d60b3329769c88388c83aa4f4579904d6a7719 doc/fix_spin_sib_npt.rst
8894d0db2c7be83e6558477cd092115c868148f76f119a1145e76502f189a433 doc/pair_spin_allegro_spin.rst
e2343cbdc0d91d81396bb8e1c6f2bca234f956ed11fce165adbf3866c60a01af doc/pair_spin_nequip_spin.rst
b8b3cceeacf32a3e0b26e5db223e1f310901244d2334a160146889e749c07187 examples/SPINGNN/README-C5.md
dba052daeff0d6b2ce61ff48a31a1c1ab42fb091221036642354dfb2fad254eb examples/SPINGNN/in.crsbr-lpsib
a91afa3276baa7956c675d756049896c659489ac4fd85b4a88aa49306fac6a8b examples/SPINGNN/in.fe-lpsib-centered
4cce5b13fb4b41e6426838bf36bb65b56762ca7c067720590ea0308d1034ac5d examples/SPINGNN/in.fe-lpsib-frozen
3b67241b833ebec72432ed7c10dff7843b1858068fae8008cea352f1d20218ed examples/SPINGNN/in.fe-lpsib-moving-nve
0beedcd7ea4207d0d60344f003455c0c31e7c863c5e06a1ff557c0f0e084f016 examples/SPINGNN/in.fe-lpsib-moving-nvt
3ac8df2322707c282d5db737881f58243397e0d6cfe31053e1979e2f6a52ec5a examples/SPINGNN/in.fe-lpsib-npt
ef62417cc80de9e17131b657fc7ad01b37304574ef8b98796550ba384bf01403 examples/SPINGNN/in.fe-sib-npt
6716ea449f24d95cc6c0b42742c48ebb02d9fd4a7cceaf42feba7d47da27ba78 examples/SPINGNN/in.fe-sib-nve
50d99a394d1e66cd268acef048a5bb8a6006659e4f21b5d65a52fd68caad74f8 examples/SPINGNN/in.fe-sib-nvt
39d5ada99edf13391eb68a3fcb39cdcd3335d87387d7ccaf91907796b440c1a1 patches/lammps.patch
aae06f610f92bece7742797fd065121741a71830869fd1f67a274c1fc6efe324 src/SPINGNN/CHANGELOG.md
bf67e41443cb5bf024d4113051f58f26f9c42ab706a2138ce3710d7fea367f3a src/SPINGNN/Install.sh
67754898a77a45ffec5e3934e6b5dcba8049b2aaf2cb215c4609af793ca10aad src/SPINGNN/KNOWN-LIMITATIONS.md
12df52ee7e77c83f6456022cc2b27e767126500d081ab85ac28bbe5843cab554 src/SPINGNN/README
0c125adbf8e9bd32062405a81074d77496ed58cfaaabafcaef94da0f8c161319 src/SPINGNN/doc-anchors.json
abaf6c22e9d5e2e900aebc4555946e67618e12cb59bf0670ca11954756d3a1e7 src/SPINGNN/fix_spin_lpsib.cpp
4d78270e1e491e5485f86646fa802bbb2849e06a7f8c932529427842a7a7f38f src/SPINGNN/fix_spin_lpsib.h
037bfb80025a8686bf4766d16ca5f2afbccdc72a84d00ffa42c3a8bcfff93add src/observe85.h
PAYLOAD_SHA256
git -C "$target_root" apply --check "$package_root/patches/lammps.patch"
mkdir -p -- "$target_root/lammps"
mkdir -p -- "$target_root/lammps/cmake"
mkdir -p -- "$target_root/lammps/cmake/Modules"
mkdir -p -- "$target_root/lammps/cmake/Modules/Packages"
mkdir -p -- "$target_root/lammps/doc"
mkdir -p -- "$target_root/lammps/doc/src"
mkdir -p -- "$target_root/lammps/examples"
mkdir -p -- "$target_root/lammps/examples/SPINGNN"
mkdir -p -- "$target_root/lammps/src"
mkdir -p -- "$target_root/lammps/src/SPIN"
mkdir -p -- "$target_root/lammps/src/SPINGNN"
mv -- "$target_root/lammps/src/SPIN/pair_allegro_spin.cpp" "$target_root/lammps/src/SPINGNN/pair_allegro_spin.cpp"
mv -- "$target_root/lammps/src/SPIN/pair_allegro_spin.h" "$target_root/lammps/src/SPINGNN/pair_allegro_spin.h"
mv -- "$target_root/lammps/src/SPIN/pair_nequip_spin.cpp" "$target_root/lammps/src/SPINGNN/pair_nequip_spin.cpp"
mv -- "$target_root/lammps/src/SPIN/pair_nequip_spin.h" "$target_root/lammps/src/SPINGNN/pair_nequip_spin.h"
cp -p -- "$package_root/cmake/SPINGNN.cmake" "$target_root/lammps/cmake/Modules/Packages/SPINGNN.cmake"
cp -p -- "$package_root/doc/SPINGNN.rst" "$target_root/lammps/doc/src/SPINGNN.rst"
cp -p -- "$package_root/doc/fix_spin_lpsib.rst" "$target_root/lammps/doc/src/fix_spin_lpsib.rst"
cp -p -- "$package_root/doc/fix_spin_sib.rst" "$target_root/lammps/doc/src/fix_spin_sib.rst"
cp -p -- "$package_root/doc/fix_spin_sib_npt.rst" "$target_root/lammps/doc/src/fix_spin_sib_npt.rst"
cp -p -- "$package_root/doc/pair_spin_allegro_spin.rst" "$target_root/lammps/doc/src/pair_spin_allegro_spin.rst"
cp -p -- "$package_root/doc/pair_spin_nequip_spin.rst" "$target_root/lammps/doc/src/pair_spin_nequip_spin.rst"
cp -p -- "$package_root/examples/SPINGNN/README-C5.md" "$target_root/lammps/examples/SPINGNN/README-C5.md"
cp -p -- "$package_root/examples/SPINGNN/in.crsbr-lpsib" "$target_root/lammps/examples/SPINGNN/in.crsbr-lpsib"
cp -p -- "$package_root/examples/SPINGNN/in.fe-lpsib-centered" "$target_root/lammps/examples/SPINGNN/in.fe-lpsib-centered"
cp -p -- "$package_root/examples/SPINGNN/in.fe-lpsib-frozen" "$target_root/lammps/examples/SPINGNN/in.fe-lpsib-frozen"
cp -p -- "$package_root/examples/SPINGNN/in.fe-lpsib-moving-nve" "$target_root/lammps/examples/SPINGNN/in.fe-lpsib-moving-nve"
cp -p -- "$package_root/examples/SPINGNN/in.fe-lpsib-moving-nvt" "$target_root/lammps/examples/SPINGNN/in.fe-lpsib-moving-nvt"
cp -p -- "$package_root/examples/SPINGNN/in.fe-lpsib-npt" "$target_root/lammps/examples/SPINGNN/in.fe-lpsib-npt"
cp -p -- "$package_root/examples/SPINGNN/in.fe-sib-npt" "$target_root/lammps/examples/SPINGNN/in.fe-sib-npt"
cp -p -- "$package_root/examples/SPINGNN/in.fe-sib-nve" "$target_root/lammps/examples/SPINGNN/in.fe-sib-nve"
cp -p -- "$package_root/examples/SPINGNN/in.fe-sib-nvt" "$target_root/lammps/examples/SPINGNN/in.fe-sib-nvt"
cp -p -- "$package_root/src/SPINGNN/CHANGELOG.md" "$target_root/lammps/src/SPINGNN/CHANGELOG.md"
cp -p -- "$package_root/src/SPINGNN/Install.sh" "$target_root/lammps/src/SPINGNN/Install.sh"
cp -p -- "$package_root/src/SPINGNN/KNOWN-LIMITATIONS.md" "$target_root/lammps/src/SPINGNN/KNOWN-LIMITATIONS.md"
cp -p -- "$package_root/src/SPINGNN/README" "$target_root/lammps/src/SPINGNN/README"
cp -p -- "$package_root/src/SPINGNN/doc-anchors.json" "$target_root/lammps/src/SPINGNN/doc-anchors.json"
cp -p -- "$package_root/src/SPINGNN/fix_spin_lpsib.cpp" "$target_root/lammps/src/SPINGNN/fix_spin_lpsib.cpp"
cp -p -- "$package_root/src/SPINGNN/fix_spin_lpsib.h" "$target_root/lammps/src/SPINGNN/fix_spin_lpsib.h"
cp -p -- "$package_root/src/observe85.h" "$target_root/lammps/src/observe85.h"
git -C "$target_root" apply "$package_root/patches/lammps.patch"
printf '%s\n' 'SPINGNN_INSTALL_OK'
