import pandas as pd

from config import (
    NEW_COTA_FILE,
    NEW_COTA_SHEET,
    OLD_ADJUDICATED_FILE,
    OLD_COTA_SHEET,
    PATIENT_ID_COLUMN,
)
from preprocessing import preprocess_new_patient
from reference_retrieval import build_reference_library, retrieve_similar_examples


def main() -> None:
    new_data = pd.read_excel(NEW_COTA_FILE, sheet_name=NEW_COTA_SHEET)
    new_data.columns = [str(column).strip() for column in new_data.columns]
    new_data["_pipeline_patient_id"] = (
        new_data[PATIENT_ID_COLUMN].replace(r"^\s*$", pd.NA, regex=True).ffill()
    )
    new_data = new_data.dropna(subset=["_pipeline_patient_id"])

    library = build_reference_library(OLD_ADJUDICATED_FILE, OLD_COTA_SHEET)
    _, first_group = next(iter(new_data.groupby("_pipeline_patient_id", sort=False)))
    patient = preprocess_new_patient(first_group.drop(columns=["_pipeline_patient_id"]))
    references = retrieve_similar_examples(patient, library, top_k=3)

    print("Setup validation passed.")
    print("Patient:", patient["patient_id"])
    print("Vendor lines:", len(patient["vendor_lines"]))
    print("Historical patients:", len(library.patients))
    print("Retrieved references:", [item["reference_patient_id"] for item in references])


if __name__ == "__main__":
    main()
