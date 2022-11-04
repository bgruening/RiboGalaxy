import sys

def main(track_type, url, sample_name, sample_description, chromosome_position, outpath):

    with open(outpath, 'w') as f:
        f. write(
f"""# copy this text into the dialog box at https://https://gwips.ucc.ie/cgi-bin/hgCustom 
# or download the file and uplaod to GWIPS on the same page
browser position {chromosome_position}
track type={track_type} name="{sample_name}" description="{sample_description}" bigDataUrl={url}""")

if __name__ == "__main__":
    track_type = sys.argv[1]
    url = sys.argv[2]
    sample_name = sys.argv[3]
    sample_description = sys.argv[4]
    chromosome_position = sys.argv[5]
    outpath = sys.argv[6]
    main(track_type, url, sample_name, sample_description, chromosome_position, outpath)
