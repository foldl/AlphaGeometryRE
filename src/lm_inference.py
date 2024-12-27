from typing import Any, Dict
from chatllm.bindings import chatllm
from chatllm.scripts import binding

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

    def __init__(self, model_file: str, mode='beam_search'):
        self.mode = mode
        self.llm = CallableLLM(chatllm.LibChatLLM(binding.PATH_BINDS), ['--hide_banner', '-m', model_file, ])

    def call(self):
        return None

    def beam_decode(
        self,
        inputs: str,
        eos_tokens: list[str]
    ):
        assert eos_tokens == [';'], "eos_tokens must be [';']"

        r = self.llm.chat(inputs).strip()

        return {
            'seqs_str': [r],
            'scores': [-1.0],
        }
