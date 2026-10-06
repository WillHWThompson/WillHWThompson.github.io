---
layout: page
permalink: /publications/
title: Publications
description: Research papers and preprints.
years: [2026,2025,2024,2023,2022,2018]
nav: true
nav_order: 1
---
<!-- _pages/publications.md -->
<div class="publications">

{%- for y in page.years %}
  <h2 class="year">{{y}}</h2>
  {% bibliography -f papers3 -q @*[year={{y}}]* %}
{% endfor %}

</div>

<div class="publications">
  <h2 id="ccm">Coherent CAPTAIN-Mills collaboration</h2>
  <p>Papers from the CCM experiment at Los Alamos National Laboratory.</p>
  {% bibliography -f ccm %}
</div>
