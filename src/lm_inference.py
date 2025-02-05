from typing import Any, Dict
from chatllm.bindings import chatllm
from chatllm.scripts import binding
import math

# Beam search parameters
BEAM_SEARCH_DEFAULT_ALPHA = 0.6
MAX_DECODE_LEN = 32

# Brevity penalty parameters
BREVITY_LEN_BIAS_NUMERATOR = 5.0
BREVITY_LEN_BIAS_DENOMINATOR = 6.0


def brevity_penalty(length: int, alpha: float = BEAM_SEARCH_DEFAULT_ALPHA):
  """Brevity penalty function for beam search penalizing short sequences.

  Args:
    alpha: float: brevity-penalty scaling parameter.
    length: int: length of considered sequence.

  Returns:
    Brevity penalty score as jax scalar.
  """
  return math.pow(
      ((BREVITY_LEN_BIAS_NUMERATOR + length) / BREVITY_LEN_BIAS_DENOMINATOR),
      alpha,
  )

class CallableLLM(chatllm.ChatLLM):

    def chat(self, user_input: str) -> str:
        self.chunk_acc = ''
        self.restart()
        super().chat(user_input)
        return self.chunk_acc

    def callback_print(self, s: str) -> None:
        self.chunk_acc = self.chunk_acc + s

class LanguageModelInference:
    """Meliad wrapper for LM inference."""

    def __init__(self, model_file: str, mode='beam_search', batch_size=2):
        self.mode = mode
        self.llm = CallableLLM(chatllm.LibChatLLM(binding.PATH_BINDS), ['--hide_banner', '-m', model_file, '--beam_size', str(batch_size)])

    def beam_decode(
        self,
        inputs: str,
        eos_tokens: list[str]
    ):
        assert eos_tokens == [';'], "eos_tokens must be [';']"

        self.llm.chat(inputs)

        results = self.llm.beam_search_results

        return {
            'seqs_str': [r['str'].strip() for r in results],
            'scores': [r['score'] / brevity_penalty(len(self.llm.text_tokenize(r['str'].strip()))) for r in results],
        }
