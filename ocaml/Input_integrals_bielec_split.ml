(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Core.Std;;

module Integrals_bielec_split : sig
(* Generate type *)
   type t = 
     {
       disk_access_ao_integrals       : Disk_access.t;
       long_range                     : bool;
       threshold_ao                   : Threshold.t;
       disk_access_mo_integrals       : Disk_access.t;
       direct                         : bool;
       disk_access_ao_integrals_standard : Disk_access.t;
       no_vvv_integrals               : bool;
       no_ivvv_integrals              : bool;
       threshold_mo                   : Threshold.t;
       mu_erf                         : float;
       disk_access_mo_integrals_standard : Disk_access.t;
       no_vvvv_integrals              : bool;
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
       disk_access_ao_integrals       : Disk_access.t;
       long_range                     : bool;
       threshold_ao                   : Threshold.t;
       disk_access_mo_integrals       : Disk_access.t;
       direct                         : bool;
       disk_access_ao_integrals_standard : Disk_access.t;
       no_vvv_integrals               : bool;
       no_ivvv_integrals              : bool;
       threshold_mo                   : Threshold.t;
       mu_erf                         : float;
       disk_access_mo_integrals_standard : Disk_access.t;
       no_vvvv_integrals              : bool;
     } with sexp
   ;;

  let get_default = Qpackage.get_ezfio_default "integrals_bielec_split";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for threshold_ao *)
  let read_threshold_ao () =
    if not (Ezfio.has_integrals_bielec_split_threshold_ao ()) then
       get_default "threshold_ao"
       |> Float.of_string
       |> Ezfio.set_integrals_bielec_split_threshold_ao
    ;
    Ezfio.get_integrals_bielec_split_threshold_ao ()
      |> Threshold.of_float
  ;;
(* Write snippet for threshold_ao *)
  let write_threshold_ao var = 
    Threshold.to_float var
    |> Ezfio.set_integrals_bielec_split_threshold_ao
  ;;

(* Read snippet for disk_access_ao_integrals *)
  let read_disk_access_ao_integrals () =
    if not (Ezfio.has_integrals_bielec_split_disk_access_ao_integrals ()) then
       get_default "disk_access_ao_integrals"
       |> String.of_string
       |> Ezfio.set_integrals_bielec_split_disk_access_ao_integrals
    ;
    Ezfio.get_integrals_bielec_split_disk_access_ao_integrals ()
      |> Disk_access.of_string
  ;;
(* Write snippet for disk_access_ao_integrals *)
  let write_disk_access_ao_integrals var = 
    Disk_access.to_string var
    |> Ezfio.set_integrals_bielec_split_disk_access_ao_integrals
  ;;

(* Read snippet for disk_access_ao_integrals_standard *)
  let read_disk_access_ao_integrals_standard () =
    if not (Ezfio.has_integrals_bielec_split_disk_access_ao_integrals_standard ()) then
       get_default "disk_access_ao_integrals_standard"
       |> String.of_string
       |> Ezfio.set_integrals_bielec_split_disk_access_ao_integrals_standard
    ;
    Ezfio.get_integrals_bielec_split_disk_access_ao_integrals_standard ()
      |> Disk_access.of_string
  ;;
(* Write snippet for disk_access_ao_integrals_standard *)
  let write_disk_access_ao_integrals_standard var = 
    Disk_access.to_string var
    |> Ezfio.set_integrals_bielec_split_disk_access_ao_integrals_standard
  ;;

(* Read snippet for disk_access_mo_integrals *)
  let read_disk_access_mo_integrals () =
    if not (Ezfio.has_integrals_bielec_split_disk_access_mo_integrals ()) then
       get_default "disk_access_mo_integrals"
       |> String.of_string
       |> Ezfio.set_integrals_bielec_split_disk_access_mo_integrals
    ;
    Ezfio.get_integrals_bielec_split_disk_access_mo_integrals ()
      |> Disk_access.of_string
  ;;
(* Write snippet for disk_access_mo_integrals *)
  let write_disk_access_mo_integrals var = 
    Disk_access.to_string var
    |> Ezfio.set_integrals_bielec_split_disk_access_mo_integrals
  ;;

(* Read snippet for disk_access_mo_integrals_standard *)
  let read_disk_access_mo_integrals_standard () =
    if not (Ezfio.has_integrals_bielec_split_disk_access_mo_integrals_standard ()) then
       get_default "disk_access_mo_integrals_standard"
       |> String.of_string
       |> Ezfio.set_integrals_bielec_split_disk_access_mo_integrals_standard
    ;
    Ezfio.get_integrals_bielec_split_disk_access_mo_integrals_standard ()
      |> Disk_access.of_string
  ;;
(* Write snippet for disk_access_mo_integrals_standard *)
  let write_disk_access_mo_integrals_standard var = 
    Disk_access.to_string var
    |> Ezfio.set_integrals_bielec_split_disk_access_mo_integrals_standard
  ;;

(* Read snippet for direct *)
  let read_direct () =
    if not (Ezfio.has_integrals_bielec_split_direct ()) then
       get_default "direct"
       |> Bool.of_string
       |> Ezfio.set_integrals_bielec_split_direct
    ;
    Ezfio.get_integrals_bielec_split_direct ()
  ;;
(* Write snippet for direct *)
  let write_direct =
     Ezfio.set_integrals_bielec_split_direct
  ;;

(* Read snippet for long_range *)
  let read_long_range () =
    if not (Ezfio.has_integrals_bielec_split_long_range ()) then
       get_default "long_range"
       |> Bool.of_string
       |> Ezfio.set_integrals_bielec_split_long_range
    ;
    Ezfio.get_integrals_bielec_split_long_range ()
  ;;
(* Write snippet for long_range *)
  let write_long_range =
     Ezfio.set_integrals_bielec_split_long_range
  ;;

(* Read snippet for threshold_mo *)
  let read_threshold_mo () =
    if not (Ezfio.has_integrals_bielec_split_threshold_mo ()) then
       get_default "threshold_mo"
       |> Float.of_string
       |> Ezfio.set_integrals_bielec_split_threshold_mo
    ;
    Ezfio.get_integrals_bielec_split_threshold_mo ()
      |> Threshold.of_float
  ;;
(* Write snippet for threshold_mo *)
  let write_threshold_mo var = 
    Threshold.to_float var
    |> Ezfio.set_integrals_bielec_split_threshold_mo
  ;;

(* Read snippet for mu_erf *)
  let read_mu_erf () =
    if not (Ezfio.has_integrals_bielec_split_mu_erf ()) then
       get_default "mu_erf"
       |> Float.of_string
       |> Ezfio.set_integrals_bielec_split_mu_erf
    ;
    Ezfio.get_integrals_bielec_split_mu_erf ()
  ;;
(* Write snippet for mu_erf *)
  let write_mu_erf =
     Ezfio.set_integrals_bielec_split_mu_erf
  ;;

(* Read snippet for no_ivvv_integrals *)
  let read_no_ivvv_integrals () =
    if not (Ezfio.has_integrals_bielec_split_no_ivvv_integrals ()) then
       get_default "no_ivvv_integrals"
       |> Bool.of_string
       |> Ezfio.set_integrals_bielec_split_no_ivvv_integrals
    ;
    Ezfio.get_integrals_bielec_split_no_ivvv_integrals ()
  ;;
(* Write snippet for no_ivvv_integrals *)
  let write_no_ivvv_integrals =
     Ezfio.set_integrals_bielec_split_no_ivvv_integrals
  ;;

(* Read snippet for no_vvv_integrals *)
  let read_no_vvv_integrals () =
    if not (Ezfio.has_integrals_bielec_split_no_vvv_integrals ()) then
       get_default "no_vvv_integrals"
       |> Bool.of_string
       |> Ezfio.set_integrals_bielec_split_no_vvv_integrals
    ;
    Ezfio.get_integrals_bielec_split_no_vvv_integrals ()
  ;;
(* Write snippet for no_vvv_integrals *)
  let write_no_vvv_integrals =
     Ezfio.set_integrals_bielec_split_no_vvv_integrals
  ;;

(* Read snippet for no_vvvv_integrals *)
  let read_no_vvvv_integrals () =
    if not (Ezfio.has_integrals_bielec_split_no_vvvv_integrals ()) then
       get_default "no_vvvv_integrals"
       |> Bool.of_string
       |> Ezfio.set_integrals_bielec_split_no_vvvv_integrals
    ;
    Ezfio.get_integrals_bielec_split_no_vvvv_integrals ()
  ;;
(* Write snippet for no_vvvv_integrals *)
  let write_no_vvvv_integrals =
     Ezfio.set_integrals_bielec_split_no_vvvv_integrals
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       disk_access_ao_integrals       = read_disk_access_ao_integrals ();
       long_range                     = read_long_range ();
       threshold_ao                   = read_threshold_ao ();
       disk_access_mo_integrals       = read_disk_access_mo_integrals ();
       direct                         = read_direct ();
       disk_access_ao_integrals_standard = read_disk_access_ao_integrals_standard ();
       no_vvv_integrals               = read_no_vvv_integrals ();
       no_ivvv_integrals              = read_no_ivvv_integrals ();
       threshold_mo                   = read_threshold_mo ();
       mu_erf                         = read_mu_erf ();
       disk_access_mo_integrals_standard = read_disk_access_mo_integrals_standard ();
       no_vvvv_integrals              = read_no_vvvv_integrals ();
     }
   ;;
(* Write all *)
   let write{ 
              disk_access_ao_integrals;
              long_range;
              threshold_ao;
              disk_access_mo_integrals;
              direct;
              disk_access_ao_integrals_standard;
              no_vvv_integrals;
              no_ivvv_integrals;
              threshold_mo;
              mu_erf;
              disk_access_mo_integrals_standard;
              no_vvvv_integrals;
            } =
     write_disk_access_ao_integrals       disk_access_ao_integrals;
     write_long_range                     long_range;
     write_threshold_ao                   threshold_ao;
     write_disk_access_mo_integrals       disk_access_mo_integrals;
     write_direct                         direct;
     write_disk_access_ao_integrals_standard disk_access_ao_integrals_standard;
     write_no_vvv_integrals               no_vvv_integrals;
     write_no_ivvv_integrals              no_ivvv_integrals;
     write_threshold_mo                   threshold_mo;
     write_mu_erf                         mu_erf;
     write_disk_access_mo_integrals_standard disk_access_mo_integrals_standard;
     write_no_vvvv_integrals              no_vvvv_integrals;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   disk_access_ao_integrals = %s
   long_range = %s
   threshold_ao = %s
   disk_access_mo_integrals = %s
   direct = %s
   disk_access_ao_integrals_standard = %s
   no_vvv_integrals = %s
   no_ivvv_integrals = %s
   threshold_mo = %s
   mu_erf = %s
   disk_access_mo_integrals_standard = %s
   no_vvvv_integrals = %s
   "
       (Disk_access.to_string b.disk_access_ao_integrals)
       (Bool.to_string b.long_range)
       (Threshold.to_string b.threshold_ao)
       (Disk_access.to_string b.disk_access_mo_integrals)
       (Bool.to_string b.direct)
       (Disk_access.to_string b.disk_access_ao_integrals_standard)
       (Bool.to_string b.no_vvv_integrals)
       (Bool.to_string b.no_ivvv_integrals)
       (Threshold.to_string b.threshold_mo)
       (Float.to_string b.mu_erf)
       (Disk_access.to_string b.disk_access_mo_integrals_standard)
       (Bool.to_string b.no_vvvv_integrals)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   Read/Write AO integrals from/to disk [ Write | Read | None ] ::
   
     disk_access_ao_integrals = %s
   
   if true, compute all the integrals using the long range interaction ::
   
     long_range = %s
   
   If |<pq|rs>| < ao_integrals_threshold then <pq|rs> is zero ::
   
     threshold_ao = %s
   
   Read/Write MO integrals from/to disk [ Write | Read | None ] ::
   
     disk_access_mo_integrals = %s
   
   Compute integrals on the fly ::
   
     direct = %s
   
   Read/Write AO integrals_standard from/to disk [ Write | Read | None ] ::
   
     disk_access_ao_integrals_standard = %s
   
   Can be switched on only if  no_vvvv_integrals  is True, then do not computes the integrals having 3 virtual orbitals ::
   
     no_vvv_integrals = %s
   
   Can be switched on only if  no_vvvv_integrals  is True, then do not computes the integrals having 3 virtual index and 1 belonging to the core inactive active orbitals ::
   
     no_ivvv_integrals = %s
   
   If |<ij|kl>| < ao_integrals_threshold then <pq|rs> is zero ::
   
     threshold_mo = %s
   
   cutting of the interaction in the range separated model ::
   
     mu_erf = %s
   
   Read/Write MO integrals_standard from/to disk [ Write | Read | None ] ::
   
     disk_access_mo_integrals_standard = %s
   
   If True, computes all integrals except for the integrals having 4 virtual index ::
   
     no_vvvv_integrals = %s
   
   "
       (Disk_access.to_string b.disk_access_ao_integrals)
       (Bool.to_string b.long_range)
       (Threshold.to_string b.threshold_ao)
       (Disk_access.to_string b.disk_access_mo_integrals)
       (Bool.to_string b.direct)
       (Disk_access.to_string b.disk_access_ao_integrals_standard)
       (Bool.to_string b.no_vvv_integrals)
       (Bool.to_string b.no_ivvv_integrals)
       (Threshold.to_string b.threshold_mo)
       (Float.to_string b.mu_erf)
       (Disk_access.to_string b.disk_access_mo_integrals_standard)
       (Bool.to_string b.no_vvvv_integrals)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end