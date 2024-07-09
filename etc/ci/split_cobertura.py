import re
from pathlib import Path


def split_cobertura(xml: Path, dst: Path) -> None:
    """
    Split the cobertura xml file into one xml file per-package and write them to dst.
    """
    dst.mkdir(exist_ok=True)
    xml_text = xml.read_text()
    minify_xml = re.compile(r"[ \t\r\n]+")
    header = xml_text[0 : xml_text.index("<packages>")]
    for i, match in enumerate(
        re.finditer(r"<package.+?</package>", xml_text, re.DOTALL)
    ):
        pkg = minify_xml.sub(" ", match[0])
        (dst / ("cobertura-%06d.xml" % i)).write_text(
            f"{header}<packages>{pkg}</packages></coverage>"
        )
