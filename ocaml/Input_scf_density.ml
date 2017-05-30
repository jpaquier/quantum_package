(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Core.Std;;

module Scf_density : sig
(* Generate type *)
   type t = 
     {
       mo_guess_type                  : MO_guess.t;
       no_oa_or_av_opt                : bool;
       level_shift                    : Positive_float.t;
       n_it_scf_max                   : Strictly_positive_int.t;
       thresh_scf                     : Threshold.t;
     } with sexp
   ;;
  val read  : unit -> t option
  val write : t-> unit
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
end = struct
(* Generate type *)
   type t = 
     {
       mo_guess_type                  : MO_guess.t;
       no_oa_or_av_opt                : bool;
       level_shift                    : Positive_float.t;
       n_it_scf_max                   : Strictly_positive_int.t;
       thresh_scf                     : Threshold.t;
     } with sexp
   ;;

  let get_default = Qpackage.get_ezfio_default "scf_density";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for level_shift *)
  let read_level_shift () =
    if not (Ezfio.has_scf_density_level_shift ()) then
       get_default "level_shift"
       |> Float.of_string
       |> Ezfio.set_scf_density_level_shift
    ;
    Ezfio.get_scf_density_level_shift ()
      |> Positive_float.of_float
  ;;
(* Write snippet for level_shift *)
  let write_level_shift var = 
    Positive_float.to_float var
    |> Ezfio.set_scf_density_level_shift
  ;;

(* Read snippet for mo_guess_type *)
  let read_mo_guess_type () =
    if not (Ezfio.has_scf_density_mo_guess_type ()) then
       get_default "mo_guess_type"
       |> String.of_string
       |> Ezfio.set_scf_density_mo_guess_type
    ;
    Ezfio.get_scf_density_mo_guess_type ()
      |> MO_guess.of_string
  ;;
(* Write snippet for mo_guess_type *)
  let write_mo_guess_type var = 
    MO_guess.to_string var
    |> Ezfio.set_scf_density_mo_guess_type
  ;;

(* Read snippet for n_it_scf_max *)
  let read_n_it_scf_max () =
    if not (Ezfio.has_scf_density_n_it_scf_max ()) then
       get_default "n_it_scf_max"
       |> Int.of_string
       |> Ezfio.set_scf_density_n_it_scf_max
    ;
    Ezfio.get_scf_density_n_it_scf_max ()
      |> Strictly_positive_int.of_int
  ;;
(* Write snippet for n_it_scf_max *)
  let write_n_it_scf_max var = 
    Strictly_positive_int.to_int var
    |> Ezfio.set_scf_density_n_it_scf_max
  ;;

(* Read snippet for no_oa_or_av_opt *)
  let read_no_oa_or_av_opt () =
    if not (Ezfio.has_scf_density_no_oa_or_av_opt ()) then
       get_default "no_oa_or_av_opt"
       |> Bool.of_string
       |> Ezfio.set_scf_density_no_oa_or_av_opt
    ;
    Ezfio.get_scf_density_no_oa_or_av_opt ()
  ;;
(* Write snippet for no_oa_or_av_opt *)
  let write_no_oa_or_av_opt =
     Ezfio.set_scf_density_no_oa_or_av_opt
  ;;

(* Read snippet for thresh_scf *)
  let read_thresh_scf () =
    if not (Ezfio.has_scf_density_thresh_scf ()) then
       get_default "thresh_scf"
       |> Float.of_string
       |> Ezfio.set_scf_density_thresh_scf
    ;
    Ezfio.get_scf_density_thresh_scf ()
      |> Threshold.of_float
  ;;
(* Write snippet for thresh_scf *)
  let write_thresh_scf var = 
    Threshold.to_float var
    |> Ezfio.set_scf_density_thresh_scf
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       mo_guess_type                  = read_mo_guess_type ();
       no_oa_or_av_opt                = read_no_oa_or_av_opt ();
       level_shift                    = read_level_shift ();
       n_it_scf_max                   = read_n_it_scf_max ();
       thresh_scf                     = read_thresh_scf ();
     }
   ;;
(* Write all *)
   let write{ 
              mo_guess_type;
              no_oa_or_av_opt;
              level_shift;
              n_it_scf_max;
              thresh_scf;
            } =
     write_mo_guess_type                  mo_guess_type;
     write_no_oa_or_av_opt                no_oa_or_av_opt;
     write_level_shift                    level_shift;
     write_n_it_scf_max                   n_it_scf_max;
     write_thresh_scf                     thresh_scf;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   mo_guess_type = %s
   no_oa_or_av_opt = %s
   level_shift = %s
   n_it_scf_max = %s
   thresh_scf = %s
   "
       (MO_guess.to_string b.mo_guess_type)
       (Bool.to_string b.no_oa_or_av_opt)
       (Positive_float.to_string b.level_shift)
       (Strictly_positive_int.to_string b.n_it_scf_max)
       (Threshold.to_string b.thresh_scf)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   Initial MO guess. Can be [ Huckel | HCore ] ::
   
     mo_guess_type = %s
   
   If true, skip the (inactive+core) --> (active) and the (active) --> (virtual) orbital rotations within the SCF procedure ::
   
     no_oa_or_av_opt = %s
   
   Energy shift on the virtual MOs to improve SCF convergence ::
   
     level_shift = %s
   
   Maximum number of SCF iterations ::
   
     n_it_scf_max = %s
   
   Threshold on the convergence of the Hartree Fock energy ::
   
     thresh_scf = %s
   
   "
       (MO_guess.to_string b.mo_guess_type)
       (Bool.to_string b.no_oa_or_av_opt)
       (Positive_float.to_string b.level_shift)
       (Strictly_positive_int.to_string b.n_it_scf_max)
       (Threshold.to_string b.thresh_scf)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end