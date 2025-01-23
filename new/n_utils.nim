# Copyright 2025 github/foldl
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

func seq2tuple2*[T](args: openArray[T]): (T, T) = (args[0], args[1])
func seq2tuple3*[T](args: openArray[T]): (T, T, T) = (args[0], args[1], args[2])
func seq2tuple4*[T](args: openArray[T]): (T, T, T, T) = (args[0], args[1], args[2], args[3])
func seq2tuple5*[T](args: openArray[T]): (T, T, T, T, T) = (args[0], args[1], args[2], args[3], args[4])
func seq2tuple6*[T](args: openArray[T]): (T, T, T, T, T, T) = (args[0], args[1], args[2], args[3], args[4], args[5])
func seq2tuple7*[T](args: openArray[T]): (T, T, T, T, T, T, T) = (args[0], args[1], args[2], args[3], args[4], args[5], args[6])
func seq2tuple8*[T](args: openArray[T]): (T, T, T, T, T, T, T, T) = (args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7])

func seq2tuple*[T](args: openArray[T]): tuple =
    case len(args)
    of 2: result = seq2tuple2(args)
    of 3: result = seq2tuple3(args)
    of 4: result = seq2tuple4(args)
    of 5: result = seq2tuple5(args)
    of 6: result = seq2tuple6(args)
    of 7: result = seq2tuple7(args)
    of 8: result = seq2tuple8(args)
    else: raise newException(ValueError, "")