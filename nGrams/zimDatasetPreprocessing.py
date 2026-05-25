import libzim
import random
from bs4 import BeautifulSoup
import re

def is_good_title(title):
    return title and title[0].isalpha()

def normalize_spacing(text):
    # Collapse all whitespace (spaces, tabs, newlines) into single space
    text = re.sub(r'\s+', ' ', text)
    return text.strip()

def extract_zim_text(zim_path, limit=40):
    archive = libzim.Archive(zim_path)
    total_entries = archive.entry_count
    print("Total entries in archive:", total_entries)
    count = 0

    # Collect valid article indices
    article_indices = []
    for i in range(total_entries):
        entry = archive._get_entry_by_id(i)
        title = entry.title
        if not is_good_title(title):
            continue
        article_indices.append(i)

    print("Valid article count:", len(article_indices))

    # Random sample
    sample_indices = random.sample(article_indices, limit)
    print("sample_indices:", sample_indices)

    for i in sample_indices:
        entry = archive._get_entry_by_id(i)

        if entry.is_redirect:
            continue

        item = entry.get_item()

        if not item.mimetype.startswith("text/html"):
            continue

        try:
            # FIX HERE
            content_bytes = bytes(item.content)
            content = content_bytes.decode("utf-8", errors="ignore")
            soup = BeautifulSoup(content, "html.parser")
            clean_text = soup.get_text(separator=" ")
            clean_preview = normalize_spacing(clean_text)

            print(f"\nTitle: {entry.title}")
            print(clean_preview[:800])  # show first 800 chars
            print("\n" + "-"*60)

        except Exception as e:
            print(f"Error reading entry {entry.title}: {e}")

zim_path = "wikipedia_es_computer_nopic_2023-10.zim"

# Example usage
extract_zim_text(zim_path)