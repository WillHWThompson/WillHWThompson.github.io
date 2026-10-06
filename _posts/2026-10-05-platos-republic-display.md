---
layout: post
title: "Plato's Republic electronic display"
date: 2026-10-05 12:00:00-0400
description: A Raspberry Pi display that reveals Plato's Republic slowly over time.
tags: making raspberry-pi
categories: projects
---



This summer, two of my friends got married (congratulations). For their wedding gift I decided to make something special. Both friends were big fans of Plato's Republic. For their gift I decided to make an electronic box that displays the
content of the republic slowly over time. The idea being that this device would sit on a desk and could be glanced at. Over time, the entirety of the Republic would be displayed.

Growing up Catholic, I was always intrigued by the lectionary calendar. The Catholic missal includes readings from the old and new testament. Depending on the day the missal contains 2 or 3 readings. The homily given by the priest usually connects or explicates the content of these readings. Over the course of a two year cycle 71% of the New Testament is read and 13.5% of the Old Testament is read. I was always struck by the idea of reading a book at this snail's pace, giving an ample amount of time to consider each and every word. I wanted to adapt this idea to another text - Plato's Republic.


I designed and laser cut the box itself out of birch sheet wood. The front is inscribed with the word "ΠΟΛΙΤΕΑ, ΠΛΑΤΩΝ" or "REPUBLIC, PLATO". The back is inscribed with a message to the bride and groom, also in Attic greek.


<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/all_unpainted.jpg" alt="Unpainted wooden pieces laid out before assembly" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/unpainted.jpg" alt="Close view of the laser-cut front panel" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/varnished_peices.jpg" alt="Varnished wooden frame pieces" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
   left all frame pieces assembled and ready for varnish. middle: closeup of the faceplate, right: varnished pieces.
</div>


I then fitted the pieces together and wood glued the frame together to let it dry.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/painted_frame.jpg" alt="Assembled and painted frame" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
   The finished and constructed frame.
</div>




The box is powered by a raspberry pi 4 with an OLED screen. The display itself is a simple webpage displayed with Chromium's kiosk mode. The changing script is made using python. The text scrolls over time so every five minutes, a new line is added to the bottom of the text and the top line is removed.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/inside.jpg" alt="Raspberry Pi and display mounted inside the frame" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
   The raspberry pi 4 and OLED screen mounted inside the frame.
</div>



<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/res_publica_photos/book_4_better_colors.jpg" alt="Finished display showing text from Plato's Republic" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
  The completed display showing republic text.
</div>
