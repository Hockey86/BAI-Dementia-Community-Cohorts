import fitz  # PyMuPDF

# Load the PDFs
doc1 = fitz.open("forest_plot_survival-BAI-Intermediate.pdf")
doc2 = fitz.open("forest_plot_survival-BAI-AgeOnly.pdf")

# Assume one page each
page1 = doc1[0]
page2 = doc2[0]

# Get sizes
width = max(page1.rect.width, page2.rect.width)
total_height = page1.rect.height + page2.rect.height+15

# Create new blank page
merged_doc = fitz.open()
new_page = merged_doc.new_page(width=width, height=total_height)

# Insert both original pages
new_page.show_pdf_page(fitz.Rect(0, page2.rect.height+15, width, total_height), doc1, 0)  # bottom
new_page.show_pdf_page(fitz.Rect(0, 3, width, page2.rect.height+3), doc2, 0)             # top

# Add labels (vector text) - multiple overlaps for bold effect
for offset in [(0,0), (0.3,0), (0,0.3), (0.3,0.3)]:
    new_page.insert_text((5+offset[0], page2.rect.height+23+offset[1]), "B. Fully adjusted",
                         fontsize=15, fontfile="/System/Library/Fonts/HelveticaNeue.ttc", fill=(0, 0, 0))
    new_page.insert_text((5+offset[0], 15+offset[1]), "A. Minimally adjusted",
                         fontsize=15, fontfile="/System/Library/Fonts/HelveticaNeue.ttc", fill=(0, 0, 0))

# Save
merged_doc.save("BAI-revision/figure3-vector.pdf")
merged_doc.close()

