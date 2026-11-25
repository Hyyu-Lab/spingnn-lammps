# Copyright (2026) Yuanmai Zheng
# SPDX-License-Identifier: GPL-2.0-only
# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# SPINGNN uses the installed SPIN base pair header.

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# The package is installed only when the SPIN base header is present.
action pair_nequip_spin.cpp pair_spin.h
action pair_nequip_spin.h pair_spin.h
action pair_allegro_spin.cpp pair_spin.h
action pair_allegro_spin.h pair_spin.h

# spin/lpsib overlay installation.
action fix_spin_lpsib.cpp pair_spin.h
action fix_spin_lpsib.h pair_spin.h

