import Lake
open Lake DSL

package photodiodeProtocol where
  version := v!"0.1.0"

@[default_target] lean_lib PhotodiodeProtocol

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.23.0"

lean_lib EntropyProofs
