import sys

def main(track_type, url, sample_name, sample_description, chromosome_position, outpath):

    def quote_attribute(value):
        if any(char in value for char in "\r\n"):
            raise ValueError("Track attributes must be single-line values")
        return value.replace("\\", "\\\\").replace('"', '\\"')

    sample_name = quote_attribute(sample_name)
    sample_description = quote_attribute(sample_description)
    url = quote_attribute(url)
    if any(char in chromosome_position for char in "\r\n"):
        raise ValueError("Browser position must be a single-line value")
    with open(outpath, 'w') as f:
        f. write(
f"""# copy this text into the dialog box at https://gwips.ucc.ie/cgi-bin/hgCustom
# or download the file and upload to GWIPS on the same page
browser position {chromosome_position}
track type={track_type} name="{sample_name}" description="{sample_description}" bigDataUrl="{url}"\n""")

if __name__ == "__main__":
    track_type = sys.argv[1]
    url = sys.argv[2]
    sample_name = sys.argv[3]
    sample_description = sys.argv[4]
    chromosome_position = sys.argv[5]
    outpath = sys.argv[6]
    main(track_type, url, sample_name, sample_description, chromosome_position, outpath)
