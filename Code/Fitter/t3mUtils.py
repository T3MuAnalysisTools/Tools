# Print statements in elagent form

class t3m_printer(_prefix):
   def __init__(self):
      self.PREFIX = _prefix

   def print(self, _str):
      print PREFIX+_str
   
   def prepend(self, _str):
      return PREFIX+_str
