from __future__ import annotations

import os
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output"

NEW_COTA_FILE = DATA_DIR / "new_cota_data.xlsx"
OLD_ADJUDICATED_FILE = DATA_DIR / "LoT Adjudication Datasets.xlsx"

NEW_COTA_SHEET = "Sheet1"
OLD_COTA_SHEET = "Cota"

PATIENT_ID_COLUMN = "cpid"
VENDOR_LOT_COLUMN = "line_of_therapy_c"
REGIMEN_COLUMN = "line_of_therapy_name"
START_DATE_COLUMN = "date_start_line_of_therapy"
END_DATE_COLUMN = "date_end_line_of_therapy"
DISCONTINUE_REASON_COLUMN = "discontinue_reason"
REFERENCE_LABEL_COLUMN = "Alberto LOT"
REFERENCE_RESIDUAL_COLUMN = "residual"

MODEL = os.getenv("OPENAI_MODEL", "gpt-5-mini")
NUMBER_OF_REFERENCE_EXAMPLES = int(os.getenv("NUMBER_OF_REFERENCE_EXAMPLES", "5"))
MAX_CONCURRENT_PATIENTS = int(os.getenv("MAX_CONCURRENT_PATIENTS", "2"))
LOW_CONFIDENCE_THRESHOLD = float(os.getenv("LOW_CONFIDENCE_THRESHOLD", "0.70"))
HIGH_CONFIDENCE_THRESHOLD = float(os.getenv("HIGH_CONFIDENCE_THRESHOLD", "0.85"))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
