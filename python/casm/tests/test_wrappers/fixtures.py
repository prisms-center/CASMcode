import os, json
from os.path import join
import casm

projects_dir = os.environ["CASM_TEST_PROJECTS"]

def read_qe_test_cases(path):
  """
  Read 'qe_test_cases.json'
  """
  return json.load(open(join(path,"qe_test_cases.json"),'r'))

def read_vasp_test_cases(path):
  """
  Read 'vasp_test_cases.json'
  """
  return json.load(open(join(path,"vasp_test_cases.json"),'r'))

def read_test_cases_eci(path):
  """
  Read 'test_cases_eci.json'
  """
  return json.load(open(join(path,"test_cases_eci.json"),'r'))