(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Core;;

module Nuclei : sig
(* Generate type *)
   type t = 
     {
       disk_access_nuclear_repulsion  : Disk_access.t;
     } [@@deriving sexp]
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
       disk_access_nuclear_repulsion  : Disk_access.t;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "nuclei";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for disk_access_nuclear_repulsion *)
  let read_disk_access_nuclear_repulsion () =
    if not (Ezfio.has_nuclei_disk_access_nuclear_repulsion ()) then
       get_default "disk_access_nuclear_repulsion"
       |> String.of_string
       |> Ezfio.set_nuclei_disk_access_nuclear_repulsion
    ;
    Ezfio.get_nuclei_disk_access_nuclear_repulsion ()
      |> Disk_access.of_string
  ;;
(* Write snippet for disk_access_nuclear_repulsion *)
  let write_disk_access_nuclear_repulsion var = 
    Disk_access.to_string var
    |> Ezfio.set_nuclei_disk_access_nuclear_repulsion
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       disk_access_nuclear_repulsion  = read_disk_access_nuclear_repulsion ();
     }
   ;;
(* Write all *)
   let write{ 
              disk_access_nuclear_repulsion;
            } =
     write_disk_access_nuclear_repulsion  disk_access_nuclear_repulsion;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   disk_access_nuclear_repulsion = %s
   "
       (Disk_access.to_string b.disk_access_nuclear_repulsion)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   Read/Write Nuclear Repulsion from/to disk [ Write | Read | None ] ::
   
     disk_access_nuclear_repulsion = %s
   
   "
       (Disk_access.to_string b.disk_access_nuclear_repulsion)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end