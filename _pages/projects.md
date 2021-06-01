---
title: "Projects"
layout: splash
permalink: /projects/
html:
  scroll-behavior: smooth
header:
  overlay_color: "#000"
  overlay_filter: "0.5"
  overlay_image: /assets/images/projects/background.jpg
excerpt: "A collection of medical tools I have built." 

intro: 
  - excerpt: <p style="font-size:14pt"> I'm passionate about building medical devices that can allow physicians to provide better and more <b> personalized care</b> for patients! The following are some of the projects I have worked on at the intersection of <b>bioengineering</b> and <b>medicine</b>! </p>

feature_row:
  - image_path: /assets/images/projects/otoai.jpg 
    alt: "OtoAI"
    title: "OtoAI"
    excerpt: "Novel digital otoscope that enables primary care physicians to take images of the inner ear and leverages machine learning to diagnose abnormal ear pathologies."
    url: "https://www.youtube.com/watch?v=oMM2NNmOSIk&t=2s"
    btn_label: "Learn More"
    btn_class: "btn--primary"
  - image_path: /assets/images/projects/kinase.jpeg
    alt: "cancer-ai"
    title: "Cancer-AI"
    excerpt: "Computational platform leveraging molecular simulation and machine learning for in silico profiling of activating kinase mutations in cancer."
    url: "https://github.com/kksuresh25/Cancer-AI"
    btn_label: "Learn More"
    btn_class: "btn--primary"
  - image_path: /assets/images/projects/tmvr.png
    alt: "mitralcrlip"
    title: "MitralClip"
    excerpt: "Interventional cardiology device for prevention of left ventricular outflow tract obstruction (LVOTO) during transcatheter mitral valve replacement."
    url: "/projects/mitralclip"
    btn_label: "Learn More"
    btn_class: "btn--primary"

---

{% include feature_row id="intro" type="center" %}

{% include feature_row %}
