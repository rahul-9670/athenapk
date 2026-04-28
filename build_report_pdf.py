"""Build PDF version of methodology_validation_full_report.md with embedded plots.

Uses reportlab for PDF generation. Parses markdown manually to handle:
- Headers (# ## ### ####)
- Plain text paragraphs
- Bullet/numbered lists
- Tables (markdown pipe-style)
- Image references
- Code blocks (triple backtick fenced)
- Inline emphasis (**bold**, *italic*, `code`)
- Horizontal rules (---)
"""

import os
import re
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.lib.colors import HexColor, black, gray, white
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, KeepTogether, HRFlowable, ListFlowable, ListItem
)
from PIL import Image as PILImage

# Configuration
DOCS_DIR = Path('/beegfs/u/bbg6470/athenapk/docs')
SRC_MD = DOCS_DIR / 'methodology_validation_full_report.md'
OUT_PDF = DOCS_DIR / 'methodology_validation_full_report.pdf'

# Page size and margins
PAGE_WIDTH, PAGE_HEIGHT = A4
LEFT_MARGIN = 2 * cm
RIGHT_MARGIN = 2 * cm
TOP_MARGIN = 2 * cm
BOT_MARGIN = 2 * cm
USABLE_WIDTH = PAGE_WIDTH - LEFT_MARGIN - RIGHT_MARGIN

# Styles
styles = getSampleStyleSheet()
STYLE_BODY = ParagraphStyle(
    name='Body', parent=styles['BodyText'],
    fontName='Helvetica', fontSize=10, leading=13,
    spaceAfter=6, alignment=TA_JUSTIFY,
)
STYLE_H1 = ParagraphStyle(
    name='H1', parent=styles['Heading1'],
    fontName='Helvetica-Bold', fontSize=18, leading=22,
    spaceAfter=12, spaceBefore=18,
    textColor=HexColor('#1a3d6e'),
)
STYLE_H2 = ParagraphStyle(
    name='H2', parent=styles['Heading2'],
    fontName='Helvetica-Bold', fontSize=14, leading=18,
    spaceAfter=8, spaceBefore=14,
    textColor=HexColor('#2a5e9c'),
)
STYLE_H3 = ParagraphStyle(
    name='H3', parent=styles['Heading3'],
    fontName='Helvetica-Bold', fontSize=12, leading=15,
    spaceAfter=6, spaceBefore=10,
    textColor=HexColor('#3a7bbf'),
)
STYLE_H4 = ParagraphStyle(
    name='H4', parent=styles['Heading4'],
    fontName='Helvetica-Bold', fontSize=11, leading=14,
    spaceAfter=4, spaceBefore=8,
    textColor=HexColor('#555555'),
)
STYLE_CODE = ParagraphStyle(
    name='Code', parent=styles['Code'],
    fontName='Courier', fontSize=8.5, leading=11,
    leftIndent=12, spaceAfter=8, spaceBefore=4,
    backColor=HexColor('#f5f5f5'),
    textColor=HexColor('#222222'),
)
STYLE_LIST = ParagraphStyle(
    name='List', parent=STYLE_BODY,
    leftIndent=18, spaceAfter=4, alignment=TA_LEFT,
)
STYLE_CAPTION = ParagraphStyle(
    name='Caption', parent=STYLE_BODY,
    fontSize=9, leading=11, alignment=TA_CENTER,
    textColor=gray, spaceAfter=12, fontName='Helvetica-Oblique',
)

def inline_format(text):
    """Convert inline markdown formatting to reportlab paragraph markup."""
    # Escape XML special chars first (but keep our markup tokens)
    text = text.replace('&', '&amp;')
    text = text.replace('<', '&lt;')
    text = text.replace('>', '&gt;')
    
    # Bold: **text**
    text = re.sub(r'\*\*([^\*]+)\*\*', r'<b>\1</b>', text)
    # Italic: *text* (but not ** which we already did)
    text = re.sub(r'(?<!\*)\*([^\*\n]+)\*(?!\*)', r'<i>\1</i>', text)
    # Inline code: `text`
    text = re.sub(r'`([^`]+)`', r'<font face="Courier" size="9" color="#553311">\1</font>', text)
    # Links [text](url) -> just keep text
    text = re.sub(r'\[([^\]]+)\]\(([^\)]+)\)', r'\1', text)
    
    return text

def make_image_flowable(img_path, max_width=USABLE_WIDTH, caption=None):
    """Build an Image flowable, scaled to fit page width."""
    img_path = str(img_path)
    if not os.path.exists(img_path):
        return Paragraph(f'<i>[missing image: {img_path}]</i>', STYLE_BODY)
    
    # Get original dimensions
    pil = PILImage.open(img_path)
    iw, ih = pil.size
    
    # Scale to fit page width (with some margin)
    target_w = min(max_width, max_width * 0.95)
    scale = target_w / iw
    target_h = ih * scale
    
    # Limit height to avoid landscape figures dominating
    max_h = (PAGE_HEIGHT - TOP_MARGIN - BOT_MARGIN) * 0.7
    if target_h > max_h:
        scale = max_h / ih
        target_h = max_h
        target_w = iw * scale
    
    img = Image(img_path, width=target_w, height=target_h)
    
    elements = [img]
    if caption:
        elements.append(Spacer(1, 4))
        elements.append(Paragraph(caption, STYLE_CAPTION))
    
    return KeepTogether(elements)

def parse_table(lines, start_idx):
    """Parse a markdown table starting at start_idx. Returns (table_flowable, end_idx)."""
    rows = []
    i = start_idx
    while i < len(lines) and lines[i].startswith('|'):
        rows.append(lines[i])
        i += 1
    
    if len(rows) < 2:
        return None, start_idx
    
    # Parse cells
    parsed = []
    for r in rows:
        # Skip separator row (---|---|...)
        if re.match(r'^\|[\s\-:|]+\|?\s*$', r):
            continue
        cells = [c.strip() for c in r.strip('|').split('|')]
        parsed.append(cells)
    
    if len(parsed) < 2:
        return None, start_idx
    
    # Build table data with formatted paragraphs
    data = []
    n_cols = len(parsed[0])
    for row_idx, row in enumerate(parsed):
        # Pad row to match column count
        while len(row) < n_cols:
            row.append('')
        row_data = []
        for cell in row[:n_cols]:
            style = ParagraphStyle(
                'cell', parent=STYLE_BODY, fontSize=8.5, leading=10,
                fontName='Helvetica-Bold' if row_idx == 0 else 'Helvetica',
                alignment=TA_LEFT,
            )
            row_data.append(Paragraph(inline_format(cell), style))
        data.append(row_data)
    
    # Compute column widths
    col_w = USABLE_WIDTH / n_cols
    
    table = Table(data, colWidths=[col_w] * n_cols, repeatRows=1)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), HexColor('#e8eef5')),
        ('GRID', (0, 0), (-1, -1), 0.4, HexColor('#888888')),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
    ]))
    return table, i

def parse_code_block(lines, start_idx):
    """Parse fenced code block. Returns (paragraph, end_idx)."""
    i = start_idx + 1
    content = []
    while i < len(lines):
        if lines[i].startswith('```'):
            i += 1
            break
        content.append(lines[i])
        i += 1
    
    code_text = '\n'.join(content)
    # Replace special chars for XML
    code_text = (code_text.replace('&', '&amp;')
                          .replace('<', '&lt;')
                          .replace('>', '&gt;')
                          .replace('\n', '<br/>')
                          .replace(' ', '&nbsp;'))
    return Paragraph(code_text, STYLE_CODE), i

def parse_markdown(md_text):
    """Parse markdown text into reportlab flowables."""
    lines = md_text.split('\n')
    flowables = []
    i = 0
    
    para_buffer = []  # accumulate paragraph lines
    
    def flush_para():
        if para_buffer:
            txt = ' '.join(para_buffer).strip()
            if txt:
                flowables.append(Paragraph(inline_format(txt), STYLE_BODY))
            para_buffer.clear()
    
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        
        # Empty line: paragraph break
        if not stripped:
            flush_para()
            i += 1
            continue
        
        # Horizontal rule
        if stripped == '---':
            flush_para()
            flowables.append(Spacer(1, 6))
            flowables.append(HRFlowable(width="100%", thickness=0.5, color=gray, spaceBefore=4, spaceAfter=8))
            i += 1
            continue
        
        # Headers
        if stripped.startswith('# '):
            flush_para()
            flowables.append(Paragraph(inline_format(stripped[2:]), STYLE_H1))
            i += 1
            continue
        if stripped.startswith('## '):
            flush_para()
            flowables.append(Paragraph(inline_format(stripped[3:]), STYLE_H2))
            i += 1
            continue
        if stripped.startswith('### '):
            flush_para()
            flowables.append(Paragraph(inline_format(stripped[4:]), STYLE_H3))
            i += 1
            continue
        if stripped.startswith('#### '):
            flush_para()
            flowables.append(Paragraph(inline_format(stripped[5:]), STYLE_H4))
            i += 1
            continue
        
        # Image: ![alt](path)
        m = re.match(r'^!\[([^\]]*)\]\(([^\)]+)\)\s*$', stripped)
        if m:
            flush_para()
            alt, path = m.groups()
            # Resolve path relative to docs/ directory
            if not os.path.isabs(path):
                path = str(DOCS_DIR / path)
            flowables.append(Spacer(1, 6))
            flowables.append(make_image_flowable(path, caption=alt if alt else None))
            flowables.append(Spacer(1, 6))
            i += 1
            continue
        
        # Code block
        if stripped.startswith('```'):
            flush_para()
            block, end_i = parse_code_block(lines, i)
            flowables.append(block)
            i = end_i
            continue
        
        # Table
        if stripped.startswith('|') and i + 1 < len(lines) and re.match(r'^\|[\s\-:|]+\|', lines[i+1].strip() if i+1 < len(lines) else ''):
            flush_para()
            tbl, end_i = parse_table(lines, i)
            if tbl:
                flowables.append(Spacer(1, 4))
                flowables.append(tbl)
                flowables.append(Spacer(1, 8))
                i = end_i
                continue
        
        # Bullet list
        if re.match(r'^[\-\*]\s+', stripped):
            flush_para()
            list_items = []
            while i < len(lines):
                m = re.match(r'^([\-\*])\s+(.*)$', lines[i].strip())
                if not m:
                    if lines[i].strip() == '':
                        # Blank line might continue list or end it - peek ahead
                        if i+1 < len(lines) and re.match(r'^[\-\*]\s+', lines[i+1].strip()):
                            i += 1
                            continue
                        else:
                            break
                    else:
                        break
                list_items.append(Paragraph(inline_format(m.group(2)), STYLE_LIST))
                i += 1
            flowables.append(ListFlowable(
                [ListItem(it, leftIndent=18) for it in list_items],
                bulletType='bullet', start='•',
                bulletFontSize=9,
            ))
            flowables.append(Spacer(1, 6))
            continue
        
        # Numbered list
        if re.match(r'^\d+\.\s+', stripped):
            flush_para()
            list_items = []
            while i < len(lines):
                m = re.match(r'^\d+\.\s+(.*)$', lines[i].strip())
                if not m:
                    if lines[i].strip() == '':
                        if i+1 < len(lines) and re.match(r'^\d+\.\s+', lines[i+1].strip()):
                            i += 1
                            continue
                        else:
                            break
                    else:
                        break
                list_items.append(Paragraph(inline_format(m.group(1)), STYLE_LIST))
                i += 1
            flowables.append(ListFlowable(
                [ListItem(it, leftIndent=18) for it in list_items],
                bulletType='1', start='1',
                bulletFontSize=9,
            ))
            flowables.append(Spacer(1, 6))
            continue
        
        # Regular paragraph line
        para_buffer.append(stripped)
        i += 1
    
    flush_para()
    return flowables

def header_footer(canvas, doc):
    """Add header and footer to each page."""
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(gray)
    # Header
    canvas.drawString(LEFT_MARGIN, PAGE_HEIGHT - 1.2*cm,
                      'AthenaPK MHD-AMR collapse methodology validation')
    canvas.drawRightString(PAGE_WIDTH - RIGHT_MARGIN, PAGE_HEIGHT - 1.2*cm,
                            'Apr 2026')
    # Footer
    canvas.drawCentredString(PAGE_WIDTH/2, 1.2*cm, f'Page {doc.page}')
    canvas.restoreState()

def main():
    print(f"Reading {SRC_MD}")
    md_text = SRC_MD.read_text(encoding='utf-8')
    print(f"  {len(md_text)} chars, {len(md_text.split(chr(10)))} lines")
    
    print("Parsing markdown...")
    flowables = parse_markdown(md_text)
    print(f"  {len(flowables)} flowables")
    
    print(f"Building PDF: {OUT_PDF}")
    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=A4,
        leftMargin=LEFT_MARGIN, rightMargin=RIGHT_MARGIN,
        topMargin=TOP_MARGIN, bottomMargin=BOT_MARGIN,
        title='AthenaPK methodology validation',
        author='Rahul Patel',
    )
    doc.build(flowables, onFirstPage=header_footer, onLaterPages=header_footer)
    
    out_size = OUT_PDF.stat().st_size
    print(f"  Done: {out_size/1024/1024:.1f} MB")

if __name__ == '__main__':
    main()
