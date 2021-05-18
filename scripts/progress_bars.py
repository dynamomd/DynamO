import abc
import math
import shutil
import sys
import time

from typing import TextIO

"""
Produce progress bar with ANSI code output.
"""
class ProgressBar(object):
    def __init__(self, target: TextIO = sys.stdout):
        self._target = target
        self._text_only = not self._target.isatty()
        self._update_width()
        self._old_progress = None
        self._old_progress_bar_str = ""

    def __enter__(self):
        if self._text_only:
            percent_str, progress_bar_str = self.generate_strings(0)
            self._target.write(progress_bar_str+"\n")
            self._target.flush()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if(exc_type is None):
            # Set to 100% for neatness, if no exception is thrown
            self.update(1.0)
        self._target.write('\n')
        self._target.flush()
        pass

    def _update_width(self):
        self._width, _ = shutil.get_terminal_size((80, 20))

    def generate_strings(self, progress : float):
        # Update width in case of resize
        self._update_width()
        # Progress bar itself
        if self._width < 12:
            # No label in excessively small terminal
            percent_str = ''
            progress_bar_str = ProgressBar.progress_bar_str(progress, self._width - 2, self._text_only)
        elif self._width < 40:
            # No padding at smaller size
            percent_str = "{:6.2f} %".format(progress * 100)
            progress_bar_str = ProgressBar.progress_bar_str(progress, self._width - 11, self._text_only) + ' '
        else:
            # Standard progress bar with padding and label
            percent_str = "{:6.2f} %".format(progress * 100) + "  "
            progress_bar_str = " " * 5 + ProgressBar.progress_bar_str(progress, self._width - 21, self._text_only)
        return percent_str, progress_bar_str
        
    def update(self, progress : float):
        if self._old_progress == progress:
            return
        self._old_progress = progress

        percent_str, progress_bar_str = self.generate_strings(progress)
        
        # Write output
        if self._text_only:
            # in text output, we just update the progress bar, missing the end parenthesis
            oldend = len(self._old_progress_bar_str[:-1].rstrip())
            newend = len(progress_bar_str[:-1].rstrip())
            #Skip the write if there's nothing to do
            if self._old_progress_bar_str == progress_bar_str:
                return
            if progress == 1:
                #We've got an update but the end is here!, write the full bar
                self._target.write(progress_bar_str[oldend:])
            else:
                self._target.write(progress_bar_str[oldend:newend])
            self._target.flush()
            self._old_progress_bar_str = progress_bar_str
        else:
            self._target.write('\033[G' + progress_bar_str + percent_str)
            self._target.flush()

    @staticmethod
    def progress_bar_str(progress : float, width : int, no_partial=False):
        # 0 <= progress <= 1
        progress = min(1, max(0, progress))
        whole_width = math.floor(progress * width)
        remainder_width = (progress * width) % 1
        part_width = math.floor(remainder_width * 8)
        if no_partial:
            part_width = 0
        part_char = [" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉"][part_width]
        if (width - whole_width - 1) < 0:
          part_char = ""
        line = "[" + "█" * whole_width + part_char + " " * (width - whole_width - 1) + "]"
        return line

