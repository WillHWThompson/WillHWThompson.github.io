---
layout: page
permalink: /publications/
title: publications
description: publications by categories in reversed chronological order.
years: [2024,2023,2022,2018]
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
