## Resubmission
This is a resubmission. In this version I have:

* Reduced the length of the title to less than 65 characters;

* Written in single quotes the BoolNet package in the DESCRIPTION;

* Replaced print() in the functions.R with message() function;

* Removed the unnecessary \dontrun{} in the example section of the 'saveDifferentiationTreeToFile' function;

* Used tempfile() and tempdir() for file saving in examples and vignettes.

## Check results
# R CMD check (local)

1 NOTE:
- checking CRAN incoming feasibility ... NOTE 
Maintainer: ‘Michele Braccini <braccini.michele@gmail.com>’ 
New submission

1 WARNING:
- ‘qpdf’ is needed for checks on size reduction of PDFs 

# devtools::check()

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

# Online checks

## devtools::check_rhub()

Platforms:
- Ubuntu Linux 20.04.1 LTS, R-release, GCC
- Fedora Linux, R-devel, clang, gfortran

checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michele Braccini <braccini.michele@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Attractor (7:32)
  Barbieri (16:38)
  Ergodic (9:29, 11:56)
  ontogenesis (13:20)
  PLOS (17:67)
  TES (2:23)
  TESs (11:40, 13:33)
  Villani (10:5, 16:27)

- Windows Server 2022, R-devel, 64 bit

checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michele Braccini <braccini.michele@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Attractor (2:23, 9:32)
  Barbieri (18:38)
  BoolNet (13:19)
  Ergodic (3:9, 11:29, 13:52)
  PLOS (19:67)
  TESs (13:36, 15:33)
  Villani (12:5, 18:27)
  ontogenesis (15:20)

checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

## devtools::check_win_devel()

1 NOTE:
- checking CRAN incoming feasibility ... NOTE 
Maintainer: 'Michele Braccini <braccini.michele@gmail.com>'
New submission

