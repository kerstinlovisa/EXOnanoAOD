import json
import os
from tabulate import tabulate

def parse_json(file_path, file_ref):
    """Parse a JSON file and extract specified fields."""
    with open(file_path, 'r', encoding='utf-8') as file:
        data = json.load(file)
    with open(file_ref, 'r', encoding='utf-8') as file:
        data_ref = json.load(file)
    groups = data['trees']["Events"]["branchgroups"]
    groups_ref = data_ref['trees']["Events"]["branchgroups"]

    nEvet = data['trees']["Events"]["entries"]
    #fields = set(groups.keys()) - set(groups_ref.keys()) 

   # Find fields that are missing in groups_ref
    new_groups = set(groups.keys()) - set(groups_ref.keys())

    # Find fields with different subs lengths (i.e. added branches in same table)
    different_subs_fields = {
        field for field in groups.keys() & groups_ref.keys()  # Only fields in both
        if len(groups[field].get("subs", [])) != len(groups_ref[field].get("subs", []))
    }
    # Combine both conditions
    fields = new_groups | different_subs_fields  

    # Extracting 'tot', 'vars', and 'items' for each field
    result = {}
    for field in fields:
        group_data = groups.get(field, {})
        group_data_ref = groups_ref.get(field, {})

        ## subtract the reference portion
        if field in different_subs_fields:
            tot =  float("%.3f"%((group_data.get("tot")-group_data_ref.get("tot"))/nEvet)) 
            var =  len(group_data.get("subs"))-len(group_data_ref.get("subs"))
        else:
            tot =  float("%.3f"%(group_data.get("tot")/nEvet))
            var =  len(group_data.get("subs"))
        result[field] = {
            "tot": tot  ,
            "vars": var ,
            "items": float("%.3f"%(group_data.get("entries")/nEvet))
        }

    result["Full nano"] = {
            "tot":float("%.3f"%(data_ref['trees']["Events"]['allsize']/nEvet)) 
    }
    
    return result

def display_results(results):
    """Format and print results as an ASCII table with multi-row headers."""
    if not results:
        print("No JSON data found.")
        return

    # Generate headers
    file_headers = ["-", "-"]  # First row with file names
    column_headers = ["Field", "vars"]  # Second row with column labels
    
    for file_name, _ in results:
        file_headers.extend([file_name, "-"] )  # Each file contributes two columns: 'tot' and 'items'
        column_headers.extend(["size(kb)/evt", "items/evt"])  # Column names per file

    # Collect all unique fields across all files
    all_fields = {field for _, data in results for field in data.keys()}

    # Prepare table data
    table = []
    total_sums = {file_name: 0 for file_name, _ in results}  # Store total sums per file
    summary_sums = {}  # Store sums for additional summary lines (e.g., DisplacedJet, etc.)

    # List of conditions for dynamic summary rows
    summary_conditions = [
        {"name": "DispJet", "pattern": "DispJet"},
        {"name": "Muon", "pattern": "Muon"},
        {"name": "MDS", "pattern": "MDS"},
        #{"name": "cscMDS", "pattern": "cscMDSCluster"},
        #{"name": "dtMDS", "pattern": "dtMDSCluster"},
    ]

    for field in sorted(all_fields):  # Sort fields alphabetically
        if "Full nano" in field: continue
        row = [field]
        vars_val = "N/A"

        # Get 'vars' value from the first occurrence of this field
        for _, data in results:
            if vars_val == "N/A" and field in data and "vars" in data[field]:
                vars_val = data[field]["vars"]

        row.append(vars_val)  # Append 'vars' (only once)

        for file_name, data in results:
            field_data = data.get(field, {})
            tot = field_data.get("tot", "N/A")
            items = field_data.get("items", "N/A")

            if tot != "N/A":
                total_sums[file_name] += float(tot)  # Sum 'tot' values per file
                #print(file_name,field,tot)
            # Check for fields matching any summary condition
                for condition in summary_conditions:
                    if condition["pattern"] in field:
                        if condition["name"] not in summary_sums:
                            summary_sums[condition["name"]] = {file_name: 0 for file_name, _ in results}
                        print(condition,file_name)
                        summary_sums[condition["name"]][file_name] += float(tot)  # Sum for matching pattern

            row.extend([tot, items])  # Append 'tot' and 'items' for each file

        table.append(row)

   # Add dynamic summary rows for each condition
    for condition in summary_conditions:
        summary_row = [condition["name"], ""]
        for file_name in summary_sums[condition["name"]]:
            summary_row.extend(["%.3f"%summary_sums[condition["name"]][file_name], ""])  # Sum only for 'tot' for the condition
        table.append(summary_row)

    # Add a "Total Sum" row at the bottom
    total_row = ["Total ExoNano sum", ""]  # No vars value for total row
    for file_name, _ in results:
        total_row.extend(["%.3f"%total_sums[file_name], ""])  # Show sum only for 'tot'
    table.append(total_row)
    
    # Add a "full nano" row at the bottom
    total_row = ["Stock nano sum", ""]  # No vars value for total row
    for _, data in results:
        total_row.extend(["%.3f"%data.get("Full nano").get("tot"), ""])  # Show sum only for 'tot'
    table.append(total_row)


    # Insert multi-row headers manually
    formatted_table = [file_headers] + [column_headers] + table

    # Print table
    print(tabulate(formatted_table, tablefmt="grid"))
 
if __name__ == "__main__":

    """Process all JSON files in a directory and extract specified fields."""
    results = []
    json_files = [  
                   #("size_EGamma.json", "size_EGamma_ref.json"),
                   ("size_EGamma_v1.json", "size_EGamma_ref.json"),
                   #("size_Muon.json","size_Muon_ref.json"),
                   ("size_Muon_v1.json","size_Muon_ref.json"),
                   #("size_TT.json","size_TT_ref.json"),
                   ("size_TT_v1.json","size_TT_ref.json")
                   ]
    for size_file, ref_file in json_files:
        parsed_data = parse_json(size_file, ref_file)
        results.append((size_file, parsed_data))

    display_results(results)
