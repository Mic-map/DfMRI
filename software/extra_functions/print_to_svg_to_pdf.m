function print_to_svg_to_pdf(figname)

print (figname, '-dsvg');
system(['rsvg-convert -f pdf -o ' figname '.pdf ' figname '.svg']);
system(['rm ' figname '.svg']);
system(['pdfcrop ' figname '.pdf ' figname '.pdf']);

end
