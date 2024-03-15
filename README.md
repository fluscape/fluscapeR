# Gravity paper 2024

## Key resources

The paper is on [overleaf](https://www.overleaf.com/project/5c3dcd9c9acf4406a2c82ae1). Please slack Steven if you need permissions in overleaf to access it.

The public code and data are part of the fluscape github organisation in the [fluscapeR R package](https://github.com/fluscape/fluscapeR).

The private code and data and the code to make the public data are at the IDDynamics github organisation in in the [fluscape repository](https://github.com/HopkinsIDD/fluscape). Please slack Justin if you need access to that.

Part of the analysis involves the creation and then use of a large binary object that contains a cache of spatial statistics used to implelent the radiation model. Recreating this rda takes a long time. A copy of the binary itself is up at Zenodo.

## Tasks

To dos
- [ ] figure out how Jon's code works to make the S matrix and rerun a new calculation that is consistent with the existing binary S matrix
- [ ] based on read_et_al_example.r write a vingette that could remake key numbers from the parameter esitmation table but would be very slow and messy
- [ ] check key numbers using that vingette for jittered and ppii versions of the data 
- [ ] put the code to actually reproduce the table into the project github, and remake the table
- [ ] remake the table fromt he jittered version of the data 
- [ ] remake one figure in the paper

Possible To Dos
- [ ] review literature for "gravity model" "radiation model" and movement ecology of SARS-CoV-2
- [ ] consider adding cites for the types of data we use in these studies
- [ ] redraft first para of intro to make sure the data make sense
- [ ] redraft conclusion

Recently completed
- [x] reproduce a radiation model result using and existing S matrix binary at the end of the gen_data.R script
- [x] reconfirm the creation of the contacts data we use for the study
- [x] organise the code on fluscape to make this easy
- [x] recreate fluscapeR into either rfluscape or fluscapeidd or immunescape [decided against]
- [x] reconfirm the jittered version of that that we add to the public database
- [x] make a really clear note of where those files are and reorganise everything else to do with gravity
- [x] make sure the reproduction is in the form of a vignette [decided against, too slow]

## Notes as of 3 October

- Submit to PLOS Comp Biol and push a preprint from there
- Check that we can recreate the key results using the read_reproduce
- Scripts will go into this directory and will call on source / R
- Will only start getting them to work after the upload

## Notes from after lunch on 3 October 2022

- Things have changed after the pandemic when lots of work was done on mobility papers

- Derek has led a review of the literature as part of the R01 they got
  - lots of work on the use of the data
  - people didn't do that much more sophisticated work
  - so this is more fundamental

## Notes on 22 November 2023

Have also found a key directory in my google drive shares. Need to just start from the very top of the results section and make sure I can remake the key outputs.

## Notes on 3 Jan 2024

Shared files are on dropbox, not google drive and this is the key file 

https://www.dropbox.com/home/shares/me_jr_hm_dc_gravity/current/to_upload?preview=read_reproduce.r

The key files in the fluscape directory are:

grav_check_and_make_upload_data.R
mob_utility_private.R
mob_utility_public.R
mob_calc_S_matrix.R
read_et_al_gen_data
read_et_al_example

## Notes on 5 Jan

Have started copying across lines from grav_check_and_make_upload_data.R to to the data-raw directory in fluscapeR. The intention is that as little code as possible for gravity is contained in the main fluscape private repo and that most of the code lives in fluscapeR. Have updated the task list here. Pretty soon, I need to dprecate this README and work from either the package of the project directory.

## Notes on 7 Jan

Still can't find the script that makes the full S matrix for fluscape, but the grav_check_and_make_upload_data.R is definitely the place to start.

## 8 Jan 2024

Searching for "s_margin" found what could likely be the right call. Now just need to make sure I can run

## 29 Jan 2023

Jon sent some code that he thinks he ran to get this to work. Need to check on that. ALso moved this file here to the fluscapeR package, rather than leaving it the manuscripts directory in the main fluscape github.

## 7 Mar 2024

SR and JR worked at FluScape lakes retreat to get the data generation working and reproduce the key model fit results. They found that the offset radiation model ran fine with the new participant and contact data from the canonical run in the fluscape repo. This was done with the old S matrix that was loaded as a binary. Jon started to look at the code to generate a new S matrix that could reproduce the same results.
