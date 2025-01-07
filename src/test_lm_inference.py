# Copyright 2023 DeepMind Technologies Limited
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

"""Unit tests for lm_inference.py."""
import os
import unittest

import test_alphageometry
import lm_inference as lm
import argparse

FLAGS = None

def parse_args():
    global FLAGS
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--model',
        type=str,
        default=':alphageometry-lm',
        help='model file name.')

    FLAGS = parser.parse_args()

class LmInferenceTest(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    super().setUpClass()

    cls.loaded_lm = lm.LanguageModelInference(
        model_file=FLAGS.model, mode='beam_search', batch_size=2
    )

  def test_lm_decode(self):
    outputs = LmInferenceTest.loaded_lm.beam_decode(
        '{S} a : ; b : ; c : ; d : T a b c d 00 T a c b d 01 ? T a d b c'
        ' {F1} x00',
        eos_tokens=[';'],
    )

    self.assertEqual(
        sorted(outputs['seqs_str']),
        sorted(['e : D a b c e 02 D a c b e 03 ;', 'e : C a c e 02 C b d e 03 ;']),
    )

if __name__ == '__main__':
  parse_args()
  test_alphageometry.run_test(LmInferenceTest)
