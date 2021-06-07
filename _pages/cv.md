---
permalink: /cv/
title: "Curriculum Vitae"
---

<!-- <object data="/assets/files/resume.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="/assets/files/resume.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="/assets/files/resume.pdf">Download PDF</a>.</p>
    </embed>
</object> -->

<!-- <object data="{{ site.url }}{{ site.baseurl }}/assets/files/resume.pdf" width="1000" height="1000" type="application/pdf"></object> -->

<!-- <iframe id="cfpdf" src="{{ site.url }}{{ site.baseurl }}/assets/files/resume.pdf" width="1000" height="1000"></iframe> -->

<iframe id="cvpdf" src="{{ site.url }}{{ site.baseurl }}/assets/files/resume.pdf" style="height:1000px; width:100%;"></iframe>
<script>
document.getElementById("theFrame").contentWindow.onload = function() {
    this.document.getElementsByTagName("img")[0].style.width="100%";
};
</script>

