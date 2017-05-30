(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Core.Std;;

module Full_ci : sig
(* Generate type *)
   type t = 
     {
       threshold_selectors_pt2        : Threshold.t;
       threshold_generators_pt2       : Threshold.t;
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
       threshold_selectors_pt2        : Threshold.t;
       threshold_generators_pt2       : Threshold.t;
     } with sexp
   ;;

  let get_default = Qpackage.get_ezfio_default "full_ci";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for threshold_generators_pt2 *)
  let read_threshold_generators_pt2 () =
    if not (Ezfio.has_full_ci_threshold_generators_pt2 ()) then
       get_default "threshold_generators_pt2"
       |> Float.of_string
       |> Ezfio.set_full_ci_threshold_generators_pt2
    ;
    Ezfio.get_full_ci_threshold_generators_pt2 ()
      |> Threshold.of_float
  ;;
(* Write snippet for threshold_generators_pt2 *)
  let write_threshold_generators_pt2 var = 
    Threshold.to_float var
    |> Ezfio.set_full_ci_threshold_generators_pt2
  ;;

(* Read snippet for threshold_selectors_pt2 *)
  let read_threshold_selectors_pt2 () =
    if not (Ezfio.has_full_ci_threshold_selectors_pt2 ()) then
       get_default "threshold_selectors_pt2"
       |> Float.of_string
       |> Ezfio.set_full_ci_threshold_selectors_pt2
    ;
    Ezfio.get_full_ci_threshold_selectors_pt2 ()
      |> Threshold.of_float
  ;;
(* Write snippet for threshold_selectors_pt2 *)
  let write_threshold_selectors_pt2 var = 
    Threshold.to_float var
    |> Ezfio.set_full_ci_threshold_selectors_pt2
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       threshold_selectors_pt2        = read_threshold_selectors_pt2 ();
       threshold_generators_pt2       = read_threshold_generators_pt2 ();
     }
   ;;
(* Write all *)
   let write{ 
              threshold_selectors_pt2;
              threshold_generators_pt2;
            } =
     write_threshold_selectors_pt2        threshold_selectors_pt2;
     write_threshold_generators_pt2       threshold_generators_pt2;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   threshold_selectors_pt2 = %s
   threshold_generators_pt2 = %s
   "
       (Threshold.to_string b.threshold_selectors_pt2)
       (Threshold.to_string b.threshold_generators_pt2)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   Thresholds on selectors (fraction of the norm) for final PT2 calculation ::
   
     threshold_selectors_pt2 = %s
   
   Thresholds on generators (fraction of the norm) for final PT2 calculation ::
   
     threshold_generators_pt2 = %s
   
   "
       (Threshold.to_string b.threshold_selectors_pt2)
       (Threshold.to_string b.threshold_generators_pt2)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end